int visEvent()
{
    // std::string hitfile = "~/alice/work/o2sim_HitsTRD_sean.root";
    std::string hitfile = "~/alice/work/sim/o2sim_HitsTRD.root";
    std::string trackletFile = "~/alice/work/sim/trdtracklets.root";
    
    TEveManager::Create();

    o2::base::GeometryManager::loadGeometry();
    auto geo = o2::trd::Geometry::instance();
    geo->createPadPlaneArray();
    geo->createClusterMatrixArray();

    int ndet = 13;
    vector<TEveBox*> vec_eve_TRD_detector_box;
    vector<TVector3> vec_TV3_local_pos;
    vec_eve_TRD_detector_box.resize(ndet);
    vec_TV3_local_pos.resize(8);

    int phosHole[19] = {402,403,404,405,406,407,432,433,434,435,436,437,462,463,464,465,466,467,538};
    int drawList[13] = {0,13,338,100,212,450,512,408,409,410,411,412,413};

    for (int idet = 0; idet < ndet; idet++)
    {
        int TRD_detector = drawList[idet];

        vec_eve_TRD_detector_box[idet] = new TEveBox;
        TString HistName = "TRD_box_";
        HistName += TRD_detector;
        vec_eve_TRD_detector_box[idet]->SetName(HistName.Data());
        vec_eve_TRD_detector_box[idet]->SetMainColor(kCyan);
        vec_eve_TRD_detector_box[idet]->SetMainTransparency(80);

        bool notInstalled = false;
        for (int missingDet : phosHole) {
            if (TRD_detector == missingDet) notInstalled = true;
        }
        if (notInstalled) continue;

        int TRD_sector = geo->getSector(TRD_detector);
        int TRD_stack = geo->getStack(TRD_detector);
        int TRD_layer = geo->getLayer(TRD_detector);
        float TRD_time0 = geo->getTime0(TRD_layer);

        float TRD_chamber_length = geo->getChamberLength(TRD_layer,TRD_stack);
        float TRD_chamber_width = geo->getChamberWidth(TRD_layer);
        float TRD_chamber_height = 8.4;

        auto padplane = geo->getPadPlane(TRD_detector);
        double TRD_row_end = padplane->getRowEnd();            // fPadRow[fNrows-1] - fLengthOPad + fPadRowSMOffset;
        double TRD_row_start = padplane->getRow0();              // fPadRow[0] + fPadRowSMOffset
        // cout << "detector: " << TRD_detector << " | " << "row start: " << TRD_row_start << endl;
        // double             TRD_row_end        = geo     ->getRowEnd(TRD_layer,TRD_stack);            // fPadRow[fNrows-1] - fLengthOPad + fPadRowSMOffset;
        // double             TRD_row_start      = geo     ->getRow0(TRD_layer,TRD_stack);              // fPadRow[0] + fPadRowSMOffset

        // Double_t             TRD_col_end        = padplane     ->GetColEnd();
        // Double_t             TRD_col_start      = padplane     ->GetCol0();
        // Double_t             TRD_row_end_ROC    = padplane     ->GetRowEndROC();         // fPadRow[fNrows-1] - fLengthOPad;
        // Double_t             TRD_col_spacing    = padplane     ->GetColSpacing();
        // Double_t             TRD_row_spacing    = padplane     ->GetRowSpacing();

        vec_TV3_local_pos[0].SetXYZ(-TRD_chamber_width/2.0, -TRD_chamber_height/2.0, -TRD_chamber_length/2.0);
        vec_TV3_local_pos[1].SetXYZ(TRD_chamber_width/2.0, -TRD_chamber_height/2.0, -TRD_chamber_length/2.0);
        vec_TV3_local_pos[2].SetXYZ(TRD_chamber_width/2.0, TRD_chamber_height/2.0, -TRD_chamber_length/2.0);
        vec_TV3_local_pos[3].SetXYZ(-TRD_chamber_width/2.0, TRD_chamber_height/2.0, -TRD_chamber_length/2.0);
        vec_TV3_local_pos[4].SetXYZ(-TRD_chamber_width/2.0, -TRD_chamber_height/2.0, TRD_chamber_length/2.0);
        vec_TV3_local_pos[5].SetXYZ(TRD_chamber_width/2.0, -TRD_chamber_height/2.0, TRD_chamber_length/2.0);
        vec_TV3_local_pos[6].SetXYZ(TRD_chamber_width/2.0, TRD_chamber_height/2.0, TRD_chamber_length/2.0);
        vec_TV3_local_pos[7].SetXYZ(-TRD_chamber_width/2.0, TRD_chamber_height/2.0, TRD_chamber_length/2.0);

        auto matrix = geo->getMatrixL2G(TRD_detector);

        for(int i_vertex = 0; i_vertex < 8; i_vertex++)
        {
            ROOT::Math::Impl::Transform3D<double>::Point loc(vec_TV3_local_pos[i_vertex][1], vec_TV3_local_pos[i_vertex][0], vec_TV3_local_pos[i_vertex][2]);
            auto glb = matrix ^ loc;

            vec_eve_TRD_detector_box[idet]->SetVertex(i_vertex, glb.x(), glb.y(), glb.z());
        }
    }

    for(int TRD_detector = 0; TRD_detector < vec_eve_TRD_detector_box.size(); TRD_detector++)
    {
        gEve->AddElement(vec_eve_TRD_detector_box[TRD_detector]);
    }
    gEve->Redraw3D();


    //-----------------------------------------------------------------------------------------------------------
    // HITS

    auto ps = new TEvePointSet();

    map<int, TEvePointSet*> psMap;
    int psColor[10] = {kRed, kBlue, kMagenta, kCyan, kOrange+7, kSpring, kTeal, kAzure, kViolet, kYellow};

    
    TFile* fin = TFile::Open(hitfile.data());
    TTree* hitTree = (TTree*)fin->Get("o2sim");

    vector<o2::trd::HitType>* hits = nullptr;
    hitTree->SetBranchAddress("TRDHit", &hits);

    int nev = hitTree->GetEntries();
    cout << "nev: " << nev << endl;;
    // nev = 1;
    
    for (int iev = 0; iev < nev; ++iev) {
        map<int, int> trackLength; 
        hitTree->GetEvent(iev);
        for (const auto& hit : *hits) {
            int det = hit.GetDetectorID();
            if (det == 408 || det == 409 || det == 410 || det == 411 || det == 412 || det == 413) {
                int trackId = hit.GetTrackID();
                trackLength[trackId]++;
            }
        }

        set<int> goodTracks;
        for (const auto& i : trackLength) {
            int length = i.second;
            if (length > 40) {
                int trackId = i.first;
                goodTracks.insert(trackId);
                psMap[trackId] = new TEvePointSet();
            }
        }

        hitTree->GetEvent(iev);
        int trackId;
        for (const auto& hit : *hits) {
            int det = hit.GetDetectorID();
            if (hit.GetTrackID() != trackId) {
                trackId = hit.GetTrackID();
            }
            bool draw = false;
            for (int drawDet=0; drawDet<ndet; drawDet++) {
                if (det == drawList[drawDet]) draw = true;
            }
            if (draw) {
                if (goodTracks.find(trackId) != goodTracks.end()) {
                    psMap[trackId]->SetNextPoint(hit.GetX(), hit.GetY(), hit.GetZ());
                    continue;
                }
                if (det == 408 || det == 409 || det == 410 || det == 411 || det == 412 || det == 413) continue;

                ps->SetNextPoint(hit.GetX(), hit.GetY(), hit.GetZ());
            }
        }
    }


    ps->SetMarkerColor(kYellow);
    ps->SetMarkerSize(1);
    ps->SetMarkerStyle(4);
    gEve->AddElement(ps);

    int colSelector = 0;
    for (const auto& x : psMap) {
        ps = x.second;
        ps->SetMarkerColor(psColor[colSelector % 10]);
        ps->SetMarkerSize(1.5);
        ps->SetMarkerStyle(4);
        gEve->AddElement(ps);
        colSelector++;
    }

    TEveLine* z_BeamLine = new TEveLine();
    z_BeamLine ->SetNextPoint(0,0,500);
    z_BeamLine ->SetNextPoint(0,0,-500);
    z_BeamLine ->SetLineColor(kMagenta);
    z_BeamLine ->SetLineWidth(2);
    gEve ->AddElement(z_BeamLine);

    auto pVertex = new TEvePointSet();
    pVertex->SetNextPoint(0,0,0);
    pVertex->SetMarkerColor(kMagenta);
    pVertex->SetMarkerSize(2.5);
    pVertex->SetMarkerStyle(4);
    gEve->AddElement(pVertex);

    gEve->Redraw3D();


    //--------------------------------------------------------------------------------------------------------------
    // TRACKLETS

    TFile* fTracklets = TFile::Open(trackletFile.data());
    TTree* trackletTree = (TTree*)fTracklets->Get("o2sim");

    vector<o2::trd::Tracklet64>* tracklets = nullptr;
    trackletTree->SetBranchAddress("Tracklet", &tracklets);

    nev = trackletTree->GetEntries();
    cout << "nev: " << nev << endl;;

    for (int iev = 0; iev < nev; ++iev) {
        trackletTree->GetEvent(iev);
        for (const auto& tracklet : *tracklets) {
            cout << tracklet.getPosition() << endl;
        }
    }

    return 0;
}