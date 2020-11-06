int vizHep()
{
    std::string hitfile = "~/alice/work/o2sim_HitsTRD.root";

    TEveManager::Create();

    auto geo = o2::trd::TRDGeometry::instance();
    geo->createPadPlaneArray();

    TGeoCombiTrans* combitrans[540];

    int ndet = 540;
    vector<TEveBox*> vec_eve_TRD_detector_box;
    vector<TVector3> vec_TV3_local_pos;
    vec_eve_TRD_detector_box.resize(ndet);
    vec_TV3_local_pos.resize(8);

    int phosHole[19] = {402,403,404,405,406,407,432,433,434,435,436,437,462,463,464,465,466,467,538};
    int drawList[13] = {0,13,338,100,212,450,512,408,409,410,411,412,413};

    for (int idet = 0; idet < ndet; idet++)
    {
        // int TRD_detector = drawList[idet];
        int TRD_detector = idet;
        vec_eve_TRD_detector_box[idet] = new TEveBox;

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
        // double             TRD_row_end        = geo     ->getRowEnd(TRD_layer,TRD_stack);            // fPadRow[fNrows-1] - fLengthOPad + fPadRowSMOffset;
        // double             TRD_row_start      = geo     ->getRow0(TRD_layer,TRD_stack);              // fPadRow[0] + fPadRowSMOffset

        // Double_t             TRD_col_end        = padplane     ->GetColEnd();
        // Double_t             TRD_col_start      = padplane     ->GetCol0();
        // Double_t             TRD_row_end_ROC    = padplane     ->GetRowEndROC();         // fPadRow[fNrows-1] - fLengthOPad;
        // Double_t             TRD_col_spacing    = padplane     ->GetColSpacing();
        // Double_t             TRD_row_spacing    = padplane     ->GetRowSpacing();

        Double_t Rotation_angle     = ((360.0/18.0)/2.0) + ((Double_t)TRD_sector)*(360.0/18.0);

        TString HistName = "TRD_box_";
        HistName += TRD_detector;
        vec_eve_TRD_detector_box[idet] ->SetName(HistName.Data());

        float sum = TRD_row_end + TRD_row_start;
        float loc[3] = {TRD_time0,0.0,float(sum/2.0)};
        float glb[3] = {0.0,0.0,0.0};
        geo ->rotateBack(TRD_detector,loc,glb);

        float locZ[3] = {float(TRD_time0-50.0),0.0,float(sum/2.0)};
        float glbZ[3] = {0.0,0.0,0.0};
        geo ->rotateBack(TRD_detector,locZ,glbZ);

        float locX[3] = {TRD_time0,50.0,float(sum/2.0)};
        float glbX[3] = {0.0,0.0,0.0};
        geo ->rotateBack(TRD_detector,locX,glbX);

        float locY[3] = {TRD_time0,0.0,float(50.0+sum/2.0)};
        float glbY[3] = {0.0,0.0,0.0};
        geo ->rotateBack(TRD_detector,locY,glbY);


        combitrans[TRD_detector] = new TGeoCombiTrans();
        combitrans[TRD_detector] ->RotateZ(Rotation_angle + 90.0);
        combitrans[TRD_detector] ->SetTranslation(glb[0],glb[1],glb[2]);

        vec_TV3_local_pos[0].SetXYZ(-TRD_chamber_width/2.0, -TRD_chamber_height/2.0, -TRD_chamber_length/2.0);
        vec_TV3_local_pos[1].SetXYZ(TRD_chamber_width/2.0, -TRD_chamber_height/2.0, -TRD_chamber_length/2.0);
        vec_TV3_local_pos[2].SetXYZ(TRD_chamber_width/2.0, TRD_chamber_height/2.0, -TRD_chamber_length/2.0);
        vec_TV3_local_pos[3].SetXYZ(-TRD_chamber_width/2.0, TRD_chamber_height/2.0, -TRD_chamber_length/2.0);
        vec_TV3_local_pos[4].SetXYZ(-TRD_chamber_width/2.0, -TRD_chamber_height/2.0, TRD_chamber_length/2.0);
        vec_TV3_local_pos[5].SetXYZ(TRD_chamber_width/2.0, -TRD_chamber_height/2.0, TRD_chamber_length/2.0);
        vec_TV3_local_pos[6].SetXYZ(TRD_chamber_width/2.0, TRD_chamber_height/2.0, TRD_chamber_length/2.0);
        vec_TV3_local_pos[7].SetXYZ(-TRD_chamber_width/2.0, TRD_chamber_height/2.0, TRD_chamber_length/2.0);

        vec_eve_TRD_detector_box[idet]->SetMainColor(kCyan);
        vec_eve_TRD_detector_box[idet]->SetMainTransparency(80);
        for(int i_vertex = 0; i_vertex < 8; i_vertex++)
        {
            Double_t arr_pos_loc[3] = {vec_TV3_local_pos[i_vertex][0], vec_TV3_local_pos[i_vertex][1], vec_TV3_local_pos[i_vertex][2]};
            Double_t arr_pos_glb[3] = {0.0,0.0,0.0};
            combitrans[TRD_detector] ->LocalToMaster(arr_pos_loc,arr_pos_glb);
            /*
            arr_pos_glb[0] += glb[0];
            arr_pos_glb[1] += glb[1];
            arr_pos_glb[2] += glb[2];
            */
            // if(TRD_detector == 106)
            // {
            //     printf("i_vertex: %d, pos: {%4.3f, %4.3f, %4.3f} \n",i_vertex,arr_pos_glb[0],arr_pos_glb[1],arr_pos_glb[2]);

            //     cout << "TRD_detector: " << TRD_detector << ", Rotation_angle: " << Rotation_angle << ", sector: " << TRD_sector
            //         << ", length: " << TRD_chamber_length << ", width: " << TRD_chamber_width << ", height: " << TRD_chamber_height
            //         << ", layer: " << TRD_layer << ", stack: " << TRD_stack << ", time0: " << TRD_time0
            //         << ", loc = {" << loc[0] << ", " << loc[1] << ", " << loc[2] << "}"
            //         << ", glb = {" << glb[0] << ", " << glb[1] << ", " << glb[2] << "}" << ", TRD_col_end: " << TRD_col_end << endl;

            // }

            vec_eve_TRD_detector_box[idet]->SetVertex(i_vertex,arr_pos_glb[0],arr_pos_glb[1],arr_pos_glb[2]);
        }
    }

    for(Int_t TRD_detector = 0; TRD_detector < vec_eve_TRD_detector_box.size(); TRD_detector++)
    {
        gEve->AddElement(vec_eve_TRD_detector_box[TRD_detector]);
    }
    gEve->Redraw3D();


    auto ps = new TEvePointSet();

    map<int, TEvePointSet*> mapPointSet;
    int psColor[11] = {kRed, kGreen, kBlue, kMagenta, kCyan, kOrange, kSpring, kTeal, kAzure, kViolet, kPink};

    
    TFile* fin = TFile::Open(hitfile.data());
    TTree* hitTree = (TTree*)fin->Get("o2sim");

    vector<o2::trd::HitType>* hits = nullptr;
    hitTree->SetBranchAddress("TRDHit", &hits);

    int nev = hitTree->GetEntries();
    // nev = 1;
    
    for (int iev = 0; iev < nev; ++iev) {
        // map<int, int> trackSizes; 
        // hitTree->GetEvent(iev);
        // for (const auto& hit : *hits) {
        //     int det = hit.GetDetectorID();
        //     if (det == 408 || det == 409 || det == 410 || det == 411 || det == 412 || det == 413) {
        //         int trackId = hit.GetTrackID();
        //         trackSizes[trackId]++;
        //     }
        // }

        // set<int> goodTracks;
        // for (const auto& i : trackSizes) {
        //     int length = i.second;
        //     if (length > 50) {
        //         int trackId = i.first;
        //         goodTracks.insert(trackId);
        //         mapPointSet[trackId] = new TEvePointSet();
        //     }
        // }

        hitTree->GetEvent(iev);
        int trackId;
        for (const auto& hit : *hits) {
            int det = hit.GetDetectorID();
            if (hit.GetTrackID() != trackId) {
                trackId = hit.GetTrackID();
            }
            bool draw = false;
            // for (int drawDet=0; drawDet<ndet; drawDet++) {
            //     if (det == drawList[drawDet]) draw = true;
            // }
            // if (draw) {
                // if (goodTracks.find(trackId) != goodTracks.end()) {
                //     mapPointSet[trackId]->SetNextPoint(hit.GetX(), hit.GetY(), hit.GetZ());
                //     continue;
                // }
                ps->SetNextPoint(hit.GetX(), hit.GetY(), hit.GetZ());

                // double locC = hit.getLocalC(); // col direction in amplification or drift volume
                // double locR = hit.getLocalR(); // row direction in amplification or drift volume
                // double locT = hit.getLocalT(); // time direction in amplification or drift volume
                // int nEl = hit.GetHitValue();
            // }
        }
    }


    ps->SetMarkerColor(kYellow);
    ps->SetMarkerSize(1);
    ps->SetMarkerStyle(4);
    gEve->AddElement(ps);

    int colSelector = 0;
    for (const auto& x : mapPointSet) {
        ps = x.second;
        ps->SetMarkerColor(psColor[colSelector % 11]);
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


    return 0;
}