/*
 * analyze.C
 * Example analysis of correlated events from (ACTAR,Mainz)-pairs of ROOT
 * files.
 */

#include <iostream>
#include <TFile.h>
#include <TTree.h>

//
// ACTAR source.
//
class ACTARSource {
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
    bool GetNext() {
      ++m_entry_i;
      if (m_entry_i == m_nentries) {
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
class MainzSource {
  public:
    MainzSource(char const *a_path, double a_slope, double a_offset):
      m_file(new TFile(a_path)),
      m_tree((TTree *)m_file->Get("triggers")),
      m_nentries(m_tree->GetEntries()),
      m_entry_i(-1),
      m_slope(a_slope),
      m_offset(a_offset),
      m_clock_time(),
      m_channel() {
        m_tree->SetBranchAddress("clock_.qtime", &m_clock_time);
        m_tree->SetBranchAddress("channel", &m_channel);
      }
    ~MainzSource() {
      delete m_file;
    }
    bool GetNext() {
      do {
        ++m_entry_i;
        if (m_entry_i == m_nentries) {
          return false;
        }
        m_tree->GetEvent(m_entry_i);
      } while (1 == m_channel);
      return true;
    }
    ULong64_t GetTS() const {
      return m_clock_time;
    }
    ULong64_t GetVulomTS() const {
      return m_slope * m_clock_time + m_offset;
    }

  private:
    TFile *m_file;
    TTree *m_tree;
    Long64_t m_nentries;
    Long64_t m_entry_i;
    double m_slope;
    double m_offset;
    ULong64_t m_clock_time;
    Int_t m_channel;
};

int main(int argc, char **argv)
{
  if (5 != argc) {
    std::cerr << "Usage: " << argv[0] << " actar.root mainz.root slope "
        "offset\n";
    exit(EXIT_FAILURE);
  }
  double slope = strtod(argv[3], NULL);
  double offset = strtod(argv[4], NULL);
  ACTARSource actar(argv[1]);
  actar.GetNext();
  MainzSource mainz(argv[2], slope, offset);
  mainz.GetNext();
  unsigned corr = 0;
  unsigned tot = 0;
  for (;;) {
    ULong64_t ts_actar = actar.GetVulomTS();
    ULong64_t ts_mainz = mainz.GetVulomTS();
    Long64_t dts = ts_actar - ts_mainz;
std::cout << ts_actar << '-' << ts_mainz << '=' << dts << '\n';
    if (abs(dts) < 2000) {
      ++corr;
      if (!actar.GetNext()) break;
      if (!mainz.GetNext()) break;
    } else if (ts_actar < ts_mainz) {
      if (!actar.GetNext()) break;
    } else {
      if (!mainz.GetNext()) break;
    }
    ++tot;
  }
  std::cout << "Correlated events = " << corr << '/' << tot << '\n';
  return 0;
}
