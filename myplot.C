{
  ifstream ifile("matcha.txt");
  if (!ifile.is_open()) {
    std::cerr << "Could not open \"matcha.txt\".\n";
    return;
  }
  int width, height;
  ifile >> width >> height;
  TH2I *th2i = new TH2I("tjo", "hej", width, 0, width, height, 0, height);
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      int score;
      ifile >> score;
      th2i->Fill(j, i, score);
    }
  }
  th2i->Draw("colz");
}
