vector<vector<double>> locMCHits; 
vector<vector<int>> locDigits; 


class SpacePointConverter
{
private:
  o2::trd::Geometry* geo;

public:
  SpacePointConverter() {
    geo = o2::trd::Geometry::instance();
    geo->createPadPlaneArray();
  }

  vector<double> calculateSpacePoint(int hcid, int padrow, int column, int position)
  {
    int idet = hcid/2;
    int ilayer = (idet % 30) % 6;
    int istack = (idet % 30) / 6;
    auto padPlane = geo->getPadPlane(ilayer, istack);

    double padWidth = padPlane->getWidthIPad();             // pad dimension in rphi-direction
    
    double padLength = padPlane->getLengthIPad();           // pad dimension in z-direction (I=inner, O=outer)
    double oPadLength = padPlane->getLengthOPad();
    
    int nPadRows = padPlane->getNrows();
   
    float driftHeight = geo->cdrHght();


    double x = driftHeight - 0.5;

    double y;
    if (hcid % 2 == 0)
    {                                       // +1.5-9 to go from center pad 10 to physical MCM center | -1 for pad 0 overlap
      y = (((double)position - 1024) * 1/75. + 1.5 - 9 - 1) * padWidth + ((double)column - 3) * padWidth * 18;     // upper edge of pad 71 at y = 0
    }
    else
    {
      y = (((double)position - 1024) * 1/75. + 1.5 + 9 - 1) * padWidth + (double)column * padWidth * 18;
    }
    
    double z = ((double)padrow - nPadRows/2 + 0.5) * padLength;          // make this the middle of the pad. this is a lower bound on z. z is between padrow*padLength and (padrow+1)*padLength
    if (padrow == nPadRows-1)
    {
      z += oPadLength - padLength;
    }
    if (padrow == 0)
    {
      z -= oPadLength - padLength;
    }

    return {x, y, z};
  }

  double calculateDy(int hcid, int slope)
  {
    int idet = hcid/2;
    int ilayer = (idet % 30) % 6;
    int istack = (idet % 30) / 6;
    auto padPlane = geo->getPadPlane(ilayer, istack);

    float driftHeight = geo->cdrHght();
    float vDrift = 1.56;
    double padWidth = padPlane->getWidthIPad();

    double dy = slope * ((driftHeight / vDrift) / 0.1) * padWidth * 1/1000.;

    return dy;
  }
};



void readMCpoints(int detector)
{
  // code borrowed from TRD/macros/checkHits.C
  // std::string hitfile = "~/alice/work/o2sim_HitsTRD_sean.root";
  std::string hitfile = "~/alice/work/sim/records/200ev/data/o2sim_HitsTRD.root";

  TFile* fin = TFile::Open(hitfile.data());
  TTree* hitTree = (TTree*)fin->Get("o2sim");

  std::vector<o2::trd::HitType>* hits = nullptr;
  hitTree->SetBranchAddress("TRDHit", &hits);

  int nev = hitTree->GetEntries();
  LOG(INFO) << nev << " MC entries found";
  // nev = 1;

  for (int iev = 0; iev < nev; ++iev) {
    hitTree->GetEntry(iev);
    for (const auto& hit : *hits) {
      int det = hit.GetDetectorID();
      if (det == detector){
        double locC = hit.getLocalC(); // col direction in amplification or drift volume
        double locR = hit.getLocalR(); // row direction in amplification or drift volume
        // LOG(info) << "locR: " << locR;
        vector<double> vtemp;
        vtemp.push_back(locC); vtemp.push_back(locR);
        locMCHits.push_back(vtemp);
      }
    }
  }
}


void readDigits(int detector)
{
  // code borrowed from TRD/macros/checkDigits.C
  // std::string digitfile = "~/alice/work/trddigits_sean.root";
  std::string digitfile = "~/alice/work/sim/records/200ev/data/trddigits.root";

  TFile* fin = TFile::Open(digitfile.data());
  TTree* digitTree = (TTree*)fin->Get("o2sim");

  std::vector<o2::trd::Digit>* digits = nullptr;
  digitTree->SetBranchAddress("TRDDigit", &digits);

  int nev = digitTree->GetEntries();
  LOG(INFO) << nev << " digit entries found";
  // nev = 1;

  for (int iev = 0; iev < nev; ++iev) {
    digitTree->GetEvent(iev);
    for (const auto& digit : *digits) {
      int det = digit.getDetector();
      if (det == detector) {
        int row = digit.getRow();
        int pad = digit.getPad();
        // LOG(info) << "row: " << row;
        vector<int> vtemp;
        vtemp.push_back(pad); vtemp.push_back(row);
        locDigits.push_back(vtemp);
      }
    }
  }
}


void plot(SpacePointConverter spc, vector<double> point, int detector)
{
  LOG(info) << "plotting";
  TCanvas* can1 = new TCanvas("trd","ptrd",0,0,1600,1600);

  gPad->Range(-75,-75,75,75);

  int ilayer = (detector % 30) % 6;
  int istack = (detector % 30) / 6;

  auto plotgeo = o2::trd::Geometry::instance();
  plotgeo->createPadPlaneArray();

  auto padPlane = plotgeo->getPadPlane(ilayer, istack);

  // inner: pads [1,142] | outer: pads 0,143
  double padWidth = padPlane->getWidthIPad();             // pad dimension in rphi-direction
  double padLength = padPlane->getLengthIPad();           // pad dimension in z-direction (I=inner, O=outer)
  double oPadWidth = padPlane->getWidthOPad();
  double oPadLength = padPlane->getLengthOPad();
  int nPadRows = padPlane->getNrows();

  // double chamberWidth = geo->getChamberWidth(ilayer);
  // double chamberLength = geo->getChamberLength(ilayer, istack);

  double padPlaneWidth = abs(padPlane->getColEnd() - padPlane->getCol0());
  double padPlaneLength = abs(padPlane->getRowEnd() - padPlane->getRow0());
  

  TBox* trdBox = new TBox(-padPlaneWidth/2, -padPlaneLength/2, padPlaneWidth/2, padPlaneLength/2);
  trdBox->SetLineWidth(1);
  trdBox->SetLineColor(kGreen);
  trdBox->Draw("l");

  TLine* vLines[144];
  TLine* hLines[nPadRows];


  // assume that padrows, columns, and MCMs each begin exactly where the previous one ends - no spaces
  // start from the middle (0,0) to take into account chamber "rim"
  // (0,0) is not the center of the physical chamber. Rather, it is the center of the pad rows and columns (alignable vol?)

  for (int i=1; i<144/2; i++){
    vLines[i] = new TLine(-padWidth*i, -padPlaneLength/2, -padWidth*i, padPlaneLength/2);
    // vLines[i]->SetLineColor(kWhite);
    if (i%18 == 0){
      vLines[i]->SetLineColor(kBlue);
    }
  }

  for (int i=144/2; i<144; i++){
    vLines[i] = new TLine(padWidth*(i-144/2), -padPlaneLength/2, padWidth*(i-144/2), padPlaneLength/2);
    // vLines[i]->SetLineColor(kWhite);
    if (i%18 == 0){
      vLines[i]->SetLineColor(kBlue);
    }
  }

  for (int i=1; i<144; i++){
    vLines[i]->Draw();
  }


  for (int i=1; i<nPadRows/2; i++){
    hLines[i] = new TLine(-padPlaneWidth/2, -padLength*i, padPlaneWidth/2, -padLength*i);
  }

  for (int i=nPadRows/2; i<nPadRows; i++){
    hLines[i] = new TLine(-padPlaneWidth/2, padLength*(i-nPadRows/2), padPlaneWidth/2, padLength*(i-nPadRows/2));
  }

  for (int i=1; i<nPadRows; i++){
    hLines[i]->Draw();
  }


  //---------------------------------------------------------------------------------------------------------------------

  TMarker* spacePoint = new TMarker(point[1], point[2], 20);
  spacePoint->SetMarkerColor(kRed);
  spacePoint->SetMarkerSize(2);
  // spacePoint->Draw();

  for (int i=0; i<locMCHits.size(); i++) {
    TMarker* mcHitMark = new TMarker(locMCHits[i][0],locMCHits[i][1],20);
    mcHitMark->SetMarkerColor(kYellow);
    mcHitMark->SetMarkerSize(1);
    mcHitMark->Draw();
  }

  for (int i=0; i<locDigits.size(); i++) {

    int pad = locDigits[i][0];
    int padrow = locDigits[i][1];
    int side = pad/72;
    int hcid = detector * 2 + side;
    int column = (pad % 72) / 18;

    int padMCM = pad % 18;
    int position = 1024 + (padMCM - 10) * 75;

    vector<double> digitPoint = spc.calculateSpacePoint(hcid, padrow, column, position);

    TMarker* digitMark = new TMarker(digitPoint[1], -digitPoint[2], 20);
    digitMark->SetMarkerColor(kRed);
    digitMark->SetMarkerSize(1);
    digitMark->Draw();
  }

  TArrow *ar1 = new TArrow(-70,-70,-70,-62,0.01,"|>");
  ar1->SetLineWidth(1.5);
  ar1->Draw();
  TArrow *ar2 = new TArrow(-70,-70,-62,-70,0.01,"|>");
  ar2->SetLineWidth(1.5);
  ar2->Draw();

  TLatex latex;
  latex.SetTextSize(0.025);
  latex.SetTextAlign(13);
  latex.DrawLatex(-68,-60,"z");
  latex.DrawLatex(-60,-68,"y");

  can1->Print("plot.pdf");
}


int hitVis()
{
  SpacePointConverter spc;

  uint64_t hcid = 677;
  uint64_t idet = hcid/2;

  // int idet = 338;
  int ilayer = (idet % 30) % 6;
  int istack = (idet % 30) / 6;

  uint64_t padrow = 3;        // 0-15
  uint64_t column = 1;        // 0-3 | always 0 at the moment
  uint64_t position = 1024;    // 0-2047 | units:pad-widths | granularity=1/75 (measured from center pad 10) 1024 is 0/center of pad 10

  vector<double> plotPoint = spc.calculateSpacePoint(hcid, padrow, column, position);


  readMCpoints(idet);

  readDigits(idet);

  plot(spc, plotPoint, idet);

  return 0;
}