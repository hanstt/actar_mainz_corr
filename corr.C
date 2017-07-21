/*
 * corr.C
 */

#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <TFile.h>
#include <TTree.h>

typedef std::pair<ULong64_t, UInt_t> TSPair;

namespace {
// # timestamps to use for startup delta-t matching.
int const c_matching_size_startup = 40;
// Used while parsing events.
int const c_matching_size_online = 20;
// Timestamp diff to accept as one event.
Long64_t const c_window = 100 * 16;
}

//
// This guy digs out noisy Heimtime timestamps.
//
class HeimtimeDigger {
  private:
    enum State {
      UNKNOWN,
      SYNC,
      ZERO_1,
      ZERO_2,
      ONE_1,
      ONE_2
    };
    static char const c_state_name[][10];
  public:
    HeimtimeDigger(double a_ms_from_ts, double a_ts_precision_ms):
      m_ms_from_ts(a_ms_from_ts),
      m_precision_ms(a_ts_precision_ms),
      m_ts_array(10),
      m_state(UNKNOWN),
      m_wip_heimtime(0),
      m_wip_i(0),
      m_my_precious() {}
    void AddHeimtimePulse(ULong64_t a_ts) {
      if (a_ts < m_ts_array.at(0)) {
        std::cerr << std::hex << "HeimtimeDigger: Out-of order timestamp! "
            "(prev=0x" << m_ts_array.at(0) << ", cur=0x" << a_ts << ")\n";
      }
      for (size_t i = m_ts_array.size() - 1; i > 0; --i) {
        m_ts_array.at(i) = m_ts_array.at(i - 1);
      }
      m_ts_array.at(0) = a_ts;

      bool has_0 = FindDT(0.16384);
      bool has_1 = FindDT(0.65536);
      bool has_sync = FindDT(5.24288);
      switch (m_state) {
        case UNKNOWN:
          if (has_sync) {
            std::cout << "HeimtimeDigger: Synced.\n";
            m_state = SYNC;
          }
          break;
        case SYNC:
          if (has_0) m_state = ZERO_1;
          else if (has_1) m_state = ONE_1;
          else if (has_sync) {
            m_wip_heimtime = 0;
            m_wip_i = 0;
          }
          else Lost();
          break;
        case ZERO_1:
          if (has_0) {
            AddBit(0, a_ts);
            m_state = ZERO_2;
          }
          else Lost();
          break;
        case ZERO_2:
          if (has_sync) m_state = SYNC;
          else Lost();
          break;
        case ONE_1:
          if (has_1) {
            AddBit(1, a_ts);
            m_state = ONE_2;
          }
          else Lost();
          break;
        case ONE_2:
          if (has_sync) m_state = SYNC;
          else Lost();
          break;
      }
    }
    // Returns timestamp and corresponding dug out Heimtime.
    TSPair GetTimestamps() const {
      return m_my_precious;
    }

  private:
    // Adds a bit to Heimtime under construction.
    void AddBit(int a_bit, ULong64_t a_ts) {
      m_wip_heimtime |= a_bit << m_wip_i;
      ++m_wip_i;
      if (32 == m_wip_i) {
        m_my_precious = std::make_pair(a_ts, m_wip_heimtime);
        m_wip_heimtime = 0;
        m_wip_i = 0;
      }
    }
    // Looks for specific duration between a pair of recent pulses.
    bool FindDT(double a_dt_ms) {
      ULong64_t t_1 = m_ts_array.at(0);
      for (size_t i = 1; i < m_ts_array.size(); ++i) {
        ULong64_t t_0 = m_ts_array.at(i);
        double dt = m_ms_from_ts * (t_1 - t_0);
        if (dt > a_dt_ms + m_precision_ms) return false;
        if (dt > a_dt_ms - m_precision_ms) return true;
      }
      return false;
    }
    void Lost() {
      std::cerr << std::dec << "HeimtimeDigger: Lost Heimtime sync (state=" <<
          c_state_name[m_state] << ").\n";
      m_state = UNKNOWN;
      m_wip_heimtime = 0;
      m_wip_i = 0;
    }

    double m_ms_from_ts;
    double m_precision_ms;
    std::vector<ULong64_t> m_ts_array;
    size_t m_ts_array_i;
    enum State m_state;
    uint32_t m_wip_heimtime;
    size_t m_wip_i;
    TSPair m_my_precious;
};
char const HeimtimeDigger::c_state_name[][10] = 
{
  "UNKNOWN",
  "SYNC",
  "ZERO_1",
  "ZERO_2",
  "ONE_1",
  "ONE_2"
};

//
// Dude syncs slope of second timestamp against first.
//
class Syncesizer {
  public:
    Syncesizer():
      m_pair_vec(50),
      m_oldest(),
      m_next(),
      m_slope(0.0) {
        Reset();
      }
    // Add pair for smooth progression.
    void Add(TSPair const &a_ts) {
      if (-1 == m_oldest) {
        m_oldest = 0;
        m_pair_vec.at(m_oldest) = a_ts;
        m_next = 1;
      } else {
        // Find slope in recent time interval.
        auto const &old = m_pair_vec.at(m_oldest);
        ULong64_t d_first = a_ts.first - old.first;
        UInt_t d_second = a_ts.second - old.second;
        m_slope = (double)d_second / d_first;
//std::cout << "Slope=" << m_slope << '\n';
        m_pair_vec.at(m_next) = a_ts;
        int next_prev = m_next;
        m_next = (m_next + 1) % m_pair_vec.size();
        if (next_prev == m_oldest) {
          m_oldest = m_next;
        }
      }
    }
    double GetSecondOverFirst() const {
      return m_slope;
    }
    void Reset() {
      m_oldest = -1;
      m_next = -1;
    }

  private:
    std::vector<TSPair> m_pair_vec;
    int m_oldest;
    int m_next;
    double m_slope;
};

//
// Delta-matches two timestamp runs.
// Pretty much Smith-Waterman, i.e. my first dynamic programming!
//
class Matcha {
  public:
    Matcha(size_t a_size, ULong64_t a_window):
      m_runs(),
      m_window(a_window) {
        for (size_t i = 0; i < 2; ++i) {
          m_runs[i].vec.resize(a_size);
        }
      }
    void Add(size_t a_i, ULong64_t a_ts) {
      assert(a_i < 2);
      auto &run = m_runs[a_i];
      run.vec.at(run.idx) = a_ts;
      run.idx = (run.idx + 1) % run.vec.size();
    }
    Long_t GetOffset() const {
      return m_offset;
    }
    bool Match() {
      // Create delta-times.
      std::vector<Long64_t> druns[2];
      for (size_t i = 0; i < 2; ++i) {
        auto const &run = m_runs[i];
        druns[i].resize(run.vec.size() - 1);
        for (size_t j = 0; j < run.vec.size() - 1; ++j) {
          size_t j0 = (j + run.idx) % run.vec.size();
          size_t j1 = (j0 + 1) % run.vec.size();
          druns[i].at(j) = run.vec.at(j1) - run.vec.at(j0);
        }
      }
      // Fill matrix and find best scoring cell.
      std::vector<Cell> matrix(druns[0].size() * druns[1].size());
      struct {
        int score;
        int i0;
        int i1;
      } best = {0, -1, -1};
      size_t ofs = 0;
      for (size_t i1 = 0; i1 < druns[1].size(); ++i1) {
        Long64_t dt1 = druns[1].at(i1);
        for (size_t i0 = 0; i0 < druns[0].size(); ++i0) {
          int const c_score_hit = 1;
          int const c_score_miss = -1;
          int const c_score_skip = -1;
          // Diagonal + pair similarity.
          Long64_t dt0 = druns[0].at(i0);
          int similarity = abs(dt1 - dt0) < m_window ? c_score_hit :
              c_score_miss;
          int prev_diag = 0 == i0 || 0 == i1 ? 0 : matrix.at(ofs - 1 -
              druns[0].size()).score;
          int diag = prev_diag + similarity;

          // Gap from left.
          int prev_left = 0 == i0 ? 0 : matrix.at(ofs - 1).score;
          int left = prev_left + c_score_skip;

          // Gap from up.
          int prev_up = 0 == i1 ? 0 : matrix.at(ofs - druns[0].size()).score;
          int up = prev_up + c_score_skip;

          // Who's the lucky winner?
          int score;
          Direction direction;
          if (diag > left) {
            if (diag > up) {
              score = diag;
              direction = DIAG;
            } else {
              score = up;
              direction = UP;
            }
          } else {
            if (left > up) {
              score = left;
              direction = LEFT;
            } else {
              score = up;
              direction = UP;
            }
          }
          score = std::max(score, 0);
          matrix.at(ofs).score = score;
          matrix.at(ofs).direction = direction;
          if (score > best.score) {
            best.score = score;
            best.i1 = i1;
            best.i0 = i0;
          }
          ++ofs;
        }
      }
      // ~50% because correlation may be split between two windows, and <50%
      // to survive spurious triggers.
      if ((unsigned)best.score < 0.4 * m_runs[0].vec.size()) {
        return false;
      }
/*      std::cout << std::dec << "Matcha: Score=" << best.score << "  (i1,i0)="
          << best.i1 << ',' << best.i0 << ".\n";*/
      ULong64_t t0 = m_runs[0].vec.at((best.i0 + m_runs[0].idx) %
          m_runs[0].vec.size());
      ULong64_t t1 = m_runs[1].vec.at((best.i1 + m_runs[1].idx) %
          m_runs[1].vec.size());
      m_offset = t0 - t1;
/*      std::cout << "Matcha: i0:" << t0 << "  i1:" << t1 << "  i0-i1:" <<
          m_offset << ".\n";
      ofstream ofile("matcha.txt");
      ofs = 0;
      ofile << druns[0].size() << ' ' << druns[1].size() << '\n';
      for (size_t i1 = 0; i1 < druns[1].size(); ++i1) {
        for (size_t i0 = 0; i0 < druns[0].size(); ++i0) {
          ofile << std::setw(4) << matrix.at(ofs).score;
          ++ofs;
        }
        ofile << '\n';
      }*/
      return true;
    }

  private:
    enum Direction {
      DIAG,
      LEFT,
      UP
    };
    struct Cell {
      int score;
      Direction direction;
    };
    struct Run {
      Run():
        vec(),
        idx(0) {}
      std::vector<ULong64_t> vec;
      size_t idx;
    };
    Run m_runs[2];
    Long64_t m_window;
    Long64_t m_offset;
};

class Source {
  public:
    virtual ~Source() {}
    virtual ULong64_t GetEntryIndex() const = 0;
    virtual bool GetNext(bool = true) = 0;
    virtual ULong64_t GetVulomTS() const = 0;
    virtual void SetEntryIndex(ULong64_t) = 0;
};

//
// ACTAR source.
//
class ACTARSource: public Source {
  public:
    ACTARSource(char const *a_path):
      m_file(new TFile(a_path)),
      m_tree((TTree *)m_file->Get("h101")),
      m_nentries(m_tree->GetEntries()),
      m_entry_i(-1),
      m_wrt(),
      m_vulom_ts() {
        m_tree->SetBranchAddress("wrt1", &m_wrt[0]);
        m_tree->SetBranchAddress("wrt2", &m_wrt[1]);
        m_tree->SetBranchAddress("wrt3", &m_wrt[2]);
        m_tree->SetBranchAddress("wrt4", &m_wrt[3]);
      }
    ~ACTARSource() {
      delete m_file;
    }
    ULong64_t GetEntryIndex() const {
      return m_entry_i;
    }
    bool GetNext(bool a_die_on_fail = true) {
      ++m_entry_i;
      if (m_entry_i == m_nentries) {
        if (a_die_on_fail) {
          std::cerr << "ACTAR file ended too fast!\n";
          exit(EXIT_FAILURE);
        }
        return false;
      }
      m_tree->GetEvent(m_entry_i);
      m_vulom_ts = (ULong64_t)m_wrt[3] << 48 | (ULong64_t)m_wrt[2] << 32 |
          (ULong64_t)m_wrt[1] << 16 | (ULong64_t)m_wrt[0];
      return true;
    }
    ULong64_t GetVulomTS() const {
      return m_vulom_ts;
    }
    void SetEntryIndex(ULong64_t a_entry_i) {
      m_entry_i = a_entry_i;
    }

  private:
    TFile *m_file;
    TTree *m_tree;
    Long64_t m_nentries;
    Long64_t m_entry_i;
    UInt_t m_wrt[4];
    ULong64_t m_vulom_ts;
};

//
// Mainz source.
//
class MainzSource: public Source {
  public:
    MainzSource(char const *a_path):
      m_file(new TFile(a_path)),
      m_tree((TTree *)m_file->Get("triggers")),
      m_nentries(m_tree->GetEntries()),
      m_entry_i(-1),
      m_clock_time(),
      m_channel(),
      m_digger(2e-6, 1e-3),
      m_syncer(),
      m_heimtime_prev(),
      m_diff0() {
        m_tree->SetBranchAddress("clock_.qtime", &m_clock_time);
        m_tree->SetBranchAddress("channel", &m_channel);
      }
    ~MainzSource() {
      delete m_file;
    }
    ULong64_t GetDiff0() const {
      return m_diff0;
    }
    ULong64_t GetEntryIndex() const {
      return m_entry_i;
    }
    ULong64_t GetHeimtime() const {
      return m_heimtime_prev;
    }
    bool GetNext(bool a_die_on_fail = true) {
      for (;;) {
        ++m_entry_i;
        if (m_entry_i == m_nentries) {
          if (a_die_on_fail) {
            std::cerr << "Main file ended too fast!\n";
            exit(EXIT_FAILURE);
          }
          return false;
        }
        m_tree->GetEvent(m_entry_i);
        if (1 == m_channel) {
          // Heimtime channel.
          m_digger.AddHeimtimePulse(m_clock_time);
          auto ts = m_digger.GetTimestamps();
          if (ts.second != m_heimtime_prev) {
            m_syncer.Add(ts);
            m_heimtime_prev = ts.second;
            if (0 == m_diff0) {
              m_diff0 = ((ULong64_t)ts.second << 28) - GetVulomTS();
            }
          }
        } else {
          break;
        }
      }
      return true;
    }
    Syncesizer const &GetSyncesizer() const {
      return m_syncer;
    }
    Syncesizer &GetSyncesizer() {
      return m_syncer;
    }
    ULong64_t GetTS() const {
      return m_clock_time;
    }
    ULong64_t GetVulomTS() const {
      return m_syncer.GetSecondOverFirst() * m_clock_time * (1 << 28) +
          m_diff0;
    }
    void SetDiff0(Long64_t a_diff0) {
      m_diff0 = a_diff0;
    }
    void SetEntryIndex(ULong64_t a_entry_i) {
      m_entry_i = a_entry_i;
    }

  private:
    TFile *m_file;
    TTree *m_tree;
    Long64_t m_nentries;
    Long64_t m_entry_i;
    ULong64_t m_clock_time;
    Int_t m_channel;
    HeimtimeDigger m_digger;
    Syncesizer m_syncer;
    UInt_t m_heimtime_prev;
    ULong64_t m_diff0;
};

int main(int argc, char **argv)
{
  if (3 != argc) {
    std::cerr << "Usage: " << argv[0] << " actar.root mainz.root\n";
    exit(EXIT_FAILURE);
  }
  ACTARSource actar(argv[1]);
  actar.GetNext();
  MainzSource mainz(argv[2]);
  mainz.GetNext();

  /*
   *
   * Perform first time synchronization.
   *
   */

  // We have 1 event from ACTAR.
  // Get 1 Heimtime from Mainz DAQ.
  while (0 == mainz.GetHeimtime()) {
    mainz.GetNext();
  }
  std::cout << std::hex << "Main: Mainz found Heimtime 0x" <<
      mainz.GetHeimtime() << ".\n";

  // Sync them up somewhat by looking for the first flip in timestamps.
  int flip = 0; // 0 = first, 1 = actar < mainz, 2 = actar > mainz.
  for (;;) {
    ULong64_t t_actar = actar.GetVulomTS();
    ULong64_t t_mainz = mainz.GetVulomTS();
    if (t_actar < t_mainz) {
      if (0 == flip) flip = 1;
      else if (2 == flip) break;
      actar.GetNext();
    } else {
      if (0 == flip) flip = 2;
      else if (1 == flip) break;
      mainz.GetNext();
    }
  }
  std::cout << std::hex << "Main: Synced ACTAR=0x" << actar.GetVulomTS() <<
      "  Mainz=0x" << mainz.GetVulomTS() << ".\n";

  // Now we should be off in the order of the alignment and transmission time
  // for a Heimtime timestamp.
  Matcha matcha(c_matching_size_startup, c_window);
  for (int k = 0; k < 2; ++k) {
    ULong64_t i_actar = actar.GetEntryIndex();
    ULong64_t i_mainz = mainz.GetEntryIndex();
    int i0, i1;
    Source *src0;
    Source *src1;
    if (0 == k) {
      i0 = k;
      i1 = 1 - k;
      src0 = &actar;
      src1 = &mainz;
    } else {
      i0 = 1 - k;
      i1 = k;
      src0 = &mainz;
      src1 = &actar;
    }
    for (size_t i = 0; i < c_matching_size_startup; ++i) {
      matcha.Add(i0, src0->GetVulomTS());
      src0->GetNext();
    }
    // Try to match delta-t against Mainz timestamp windows.
    for (;;) {
      for (size_t i = 0; i < c_matching_size_startup; ++i) {
        matcha.Add(i1, src1->GetVulomTS());
        src1->GetNext();
      }
      // Cross your fingers.
      if (matcha.Match()) {
        goto match_found;
      }
      // If we're getting way out of hand (>1s), then there's simply no sync.
      if (src1->GetVulomTS() > src0->GetVulomTS() + 16e8) {
        break;
      }
    }
    actar.SetEntryIndex(i_actar - 1);
    actar.GetNext();
    mainz.SetEntryIndex(i_mainz - 1);
    mainz.GetNext();
  }
  std::cerr << "No sync found!\n";
  exit(EXIT_FAILURE);
match_found:
  mainz.SetDiff0(mainz.GetDiff0() + matcha.GetOffset());

  /*
   *
   * Main loop with smaller and faster synchronization.
   *
   */

  actar.SetEntryIndex(-1);
  actar.GetNext();
  mainz.SetEntryIndex(-1);
  mainz.GetNext();
  mainz.GetSyncesizer().Reset();
  Matcha matcha_online(c_matching_size_online, c_window);
  unsigned hits = 0;
  for (unsigned counter = 0;;) {
    ULong64_t ts_actar = actar.GetVulomTS();
    ULong64_t ts_mainz = mainz.GetVulomTS();
    Long64_t dt = ts_actar - ts_mainz;
    dt = 0 > dt ? -dt : dt;
std::cout << ts_actar << '-' << ts_mainz << '=' << dt << '\n';
    bool do_next_actar = false;
    bool do_next_mainz = false;
    if (dt < c_window) {
      // Event!
      ++hits;
      do_next_actar = true;
      do_next_mainz = true;
    } else if (ts_actar < ts_mainz) {
      do_next_actar = true;
    } else {
      do_next_mainz = true;
    }
    if (do_next_actar) {
      if (!actar.GetNext(false)) {
        break;
      }
      matcha_online.Add(0, actar.GetVulomTS());
    }
    if (do_next_mainz) {
      if (!mainz.GetNext(false)) {
        break;
      }
      matcha_online.Add(1, mainz.GetVulomTS());
    }
    counter = (counter + 1) % c_matching_size_online;
    if (0 == counter) {
      if (matcha_online.Match()) {
        mainz.SetDiff0(mainz.GetDiff0() + matcha_online.GetOffset());
      }
    }
  }
  std::cout << "Well, that ended well!\n";
  return 0;
}
