struct tree_variables
{
    int run, lumis, event;
    float mgg, mjj, mtot;
    float costhetastar_CS;
    float jj_DR, gg_DR;
    float pho1_pt, pho1_e, pho1_phi, pho1_eta, pho1_mass;
    float pho2_pt, pho2_e, pho2_phi, pho2_eta, pho2_mass;
    float pho1_r9, pho2_r9;
    float jet1_pt, jet1_e, jet1_phi, jet1_eta, jet1_mass, jet1_btagSF_M, jet1_btagEff_M, jet1_btagSFErrorUp_M, jet1_btagSFErrorDown_M, jet1_btagEffError_M, jet1_csvBtag;
    float jet2_pt, jet2_e, jet2_phi, jet2_eta, jet2_mass, jet2_btagSF_M, jet2_btagEff_M, jet2_btagSFErrorUp_M, jet2_btagSFErrorDown_M, jet2_btagEffError_M, jet2_csvBtag;
    int jet1_flavour;
    int jet2_flavour;
    float mjj_wokinfit, mtot_wokinfit, mjj_wkinfit;
    int cut_based_ct, njets_kRadionID_and_CSVM, selection_cut_level;
    float weight, evWeight, evWeight_w_btagSF;
    float regcosthetastar, minDRgregkinj;
    int njets_kLooseID;

    float weightBtagSF, weightBtagSFerrUp, weightBtagSFerrDown;
};

void setup_intree(TTree *intree, tree_variables *t, std::string whichJet)
{
    intree->SetBranchAddress("run", &t->run);
    intree->SetBranchAddress("lumis", &t->lumis);
    intree->SetBranchAddress("event", &t->event);
    intree->SetBranchAddress("gg_mass", &t->mgg);
    intree->SetBranchAddress("pho1_pt", &t->pho1_pt);
    intree->SetBranchAddress("pho1_e", &t->pho1_e);
    intree->SetBranchAddress("pho1_phi", &t->pho1_phi);
    intree->SetBranchAddress("pho1_eta", &t->pho1_eta);
    intree->SetBranchAddress("pho1_mass", &t->pho1_mass);
    intree->SetBranchAddress("pho1_r9", &t->pho1_r9);
    intree->SetBranchAddress("pho2_pt", &t->pho2_pt);
    intree->SetBranchAddress("pho2_e", &t->pho2_e);
    intree->SetBranchAddress("pho2_phi", &t->pho2_phi);
    intree->SetBranchAddress("pho2_eta", &t->pho2_eta);
    intree->SetBranchAddress("pho2_mass", &t->pho2_mass);
    intree->SetBranchAddress("pho2_r9", &t->pho2_r9);
    intree->SetBranchAddress(Form("%sjet1_pt", whichJet.c_str()), &t->jet1_pt);
    intree->SetBranchAddress(Form("%sjet1_e", whichJet.c_str()), &t->jet1_e);
    intree->SetBranchAddress(Form("%sjet1_phi", whichJet.c_str()), &t->jet1_phi);
    intree->SetBranchAddress(Form("%sjet1_eta", whichJet.c_str()), &t->jet1_eta);
    intree->SetBranchAddress(Form("%sjet1_mass", whichJet.c_str()), &t->jet1_mass);
    intree->SetBranchAddress(Form("%sjet1_csvBtag", whichJet.c_str()), &t->jet1_csvBtag);
    intree->SetBranchAddress(Form("%sjet1_btagSF_M", whichJet.c_str()), &t->jet1_btagSF_M);
    intree->SetBranchAddress(Form("%sjet1_btagSFErrorUp_M", whichJet.c_str()), &t->jet1_btagSFErrorUp_M);
    intree->SetBranchAddress(Form("%sjet1_btagSFErrorDown_M", whichJet.c_str()), &t->jet1_btagSFErrorDown_M);
    intree->SetBranchAddress(Form("%sjet1_btagEff_M", whichJet.c_str()), &t->jet1_btagEff_M);
    intree->SetBranchAddress(Form("%sjet1_btagEffError_M", whichJet.c_str()), &t->jet1_btagEffError_M);
    intree->SetBranchAddress(Form("%sjet1_flavour", whichJet.c_str()), &t->jet1_flavour);
    intree->SetBranchAddress(Form("%sjet2_pt", whichJet.c_str()), &t->jet2_pt);
    intree->SetBranchAddress(Form("%sjet2_e", whichJet.c_str()), &t->jet2_e);
    intree->SetBranchAddress(Form("%sjet2_phi", whichJet.c_str()), &t->jet2_phi);
    intree->SetBranchAddress(Form("%sjet2_eta", whichJet.c_str()), &t->jet2_eta);
    intree->SetBranchAddress(Form("%sjet2_mass", whichJet.c_str()), &t->jet2_mass);
    intree->SetBranchAddress(Form("%sjet2_csvBtag", whichJet.c_str()), &t->jet2_csvBtag);
    intree->SetBranchAddress(Form("%sjet2_btagSF_M", whichJet.c_str()), &t->jet2_btagSF_M);
    intree->SetBranchAddress(Form("%sjet2_btagSFErrorUp_M", whichJet.c_str()), &t->jet2_btagSFErrorUp_M);
    intree->SetBranchAddress(Form("%sjet2_btagSFErrorDown_M", whichJet.c_str()), &t->jet2_btagSFErrorDown_M);
    intree->SetBranchAddress(Form("%sjet2_btagEff_M", whichJet.c_str()), &t->jet2_btagEff_M);
    intree->SetBranchAddress(Form("%sjet2_btagEffError_M", whichJet.c_str()), &t->jet2_btagEffError_M);
    intree->SetBranchAddress(Form("%sjet2_flavour", whichJet.c_str()), &t->jet2_flavour);
    intree->SetBranchAddress(Form("%sjj_mass", whichJet.c_str()), &t->mjj);
    intree->SetBranchAddress(Form("%sggjj_mass", whichJet.c_str()), &t->mtot);
    intree->SetBranchAddress(Form("%scosthetastar_CS", whichJet.c_str()), &t->costhetastar_CS);
    intree->SetBranchAddress(Form("%sjj_DR", whichJet.c_str()), &t->jj_DR);
    intree->SetBranchAddress("gg_DR", &t->gg_DR);
// Prepare mjj and mggjj variables "without kin fit" on which to cut
// (in case there is no kin fit asked for, they are just a dumb copy/paste)
    if( (strcmp("kin", whichJet.c_str()) == 0) || (strcmp("regkin", whichJet.c_str()) == 0) )
    {
        std::string whichJet_tmp = "";
        if(strcmp("regkin", whichJet.c_str()) == 0) whichJet_tmp = "reg";
        intree->SetBranchAddress(Form("%sjj_mass", whichJet_tmp.c_str()), &t->mjj_wokinfit);
        intree->SetBranchAddress(Form("%sggjj_mass", whichJet_tmp.c_str()), &t->mtot_wokinfit);
    }
    intree->SetBranchAddress("njets_kRadionID_and_CSVM", &t->njets_kRadionID_and_CSVM);
    intree->SetBranchAddress("selection_cut_level", &t->selection_cut_level);
    intree->SetBranchAddress("evweight", &t->evWeight);
    intree->SetBranchAddress("weight", &t->weight);
    intree->SetBranchAddress("regcosthetastar", &t->regcosthetastar);
    intree->SetBranchAddress("minDRgregkinj", &t->minDRgregkinj);
    intree->SetBranchAddress("njets_kLooseID", &t->njets_kLooseID);
    return;
}

void setup_outtree(TTree *outtree, tree_variables *t)
{
    outtree->Branch("run", &t->run, "run/I");
    outtree->Branch("lumis", &t->lumis, "lumis/I");
    outtree->Branch("event", &t->event, "event/I");
    outtree->Branch("pho1_pt", &t->pho1_pt, "pho1_pt/F");
    outtree->Branch("pho1_e", &t->pho1_e, "pho1_e/F");
    outtree->Branch("pho1_phi", &t->pho1_phi, "pho1_phi/F");
    outtree->Branch("pho1_eta", &t->pho1_eta, "pho1_eta/F");
    outtree->Branch("pho1_mass", &t->pho1_mass, "pho1_mass/F");
    outtree->Branch("pho1_r9", &t->pho1_r9, "pho1_r9/F");
    outtree->Branch("pho2_pt", &t->pho2_pt, "pho2_pt/F");
    outtree->Branch("pho2_e", &t->pho2_e, "pho2_e/F");
    outtree->Branch("pho2_phi", &t->pho2_phi, "pho2_phi/F");
    outtree->Branch("pho2_eta", &t->pho2_eta, "pho2_eta/F");
    outtree->Branch("pho2_mass", &t->pho2_mass, "pho2_mass/F");
    outtree->Branch("pho2_r9", &t->pho2_r9, "pho2_r9/F");
    outtree->Branch("jet1_pt", &t->jet1_pt, "jet1_pt/F");
    outtree->Branch("jet1_e", &t->jet1_e, "jet1_e/F");
    outtree->Branch("jet1_phi", &t->jet1_phi, "jet1_phi/F");
    outtree->Branch("jet1_eta", &t->jet1_eta, "jet1_eta/F");
    outtree->Branch("jet1_mass", &t->jet1_mass, "jet1_mass/F");
    outtree->Branch("jet1_csvBtag", &t->jet1_csvBtag, "jet1_csvBtag/F");
    outtree->Branch("jet1_btagSF_M", &t->jet1_btagSF_M, "jet1_btagSF_M/F");
    outtree->Branch("jet1_btagSFErrorUp_M", &t->jet1_btagSFErrorUp_M, "jet1_btagSFErrorUp_M/F");
    outtree->Branch("jet1_btagSFErrorDown_M", &t->jet1_btagSFErrorDown_M, "jet1_btagSFErrorDown_M/F");
    outtree->Branch("jet1_btagEff_M", &t->jet1_btagEff_M, "jet1_btagEff_M/F");
    outtree->Branch("jet1_btagEffError_M", &t->jet1_btagEffError_M, "jet1_btagEffError_M/F");
    outtree->Branch("jet1_flavour", &t->jet1_flavour, "jet1_flavour/I");
    outtree->Branch("jet2_pt", &t->jet2_pt, "jet2_pt/F");
    outtree->Branch("jet2_e", &t->jet2_e, "jet2_e/F");
    outtree->Branch("jet2_phi", &t->jet2_phi, "jet2_phi/F");
    outtree->Branch("jet2_eta", &t->jet2_eta, "jet2_eta/F");
    outtree->Branch("jet2_mass", &t->jet2_mass, "jet2_mass/F");
    outtree->Branch("jet2_csvBtag", &t->jet2_csvBtag, "jet2_csvBtag/F");
    outtree->Branch("jet2_btagSF_M", &t->jet2_btagSF_M, "jet2_btagSF_M/F");
    outtree->Branch("jet2_btagSFErrorUp_M", &t->jet2_btagSFErrorUp_M, "jet2_btagSFErrorUp_M/F");
    outtree->Branch("jet2_btagSFErrorDown_M", &t->jet2_btagSFErrorDown_M, "jet2_btagSFErrorDown_M/F");
    outtree->Branch("jet2_btagEff_M", &t->jet2_btagEff_M, "jet2_btagEff_M/F");
    outtree->Branch("jet2_btagEffError_M", &t->jet2_btagEffError_M, "jet2_btagEffError_M/F");
    outtree->Branch("jet2_flavour", &t->jet2_flavour, "jet2_flavour/I");
    outtree->Branch("mgg", &t->mgg, "mgg/F");
    outtree->Branch("mjj", &t->mjj, "mjj/F");
    outtree->Branch("mtot", &t->mtot, "mtot/F");
    outtree->Branch("costhetastar_CS", &t->costhetastar_CS, "costhetastar_CS/F");
    outtree->Branch("jj_DR", &t->jj_DR, "jj_DR/F");
    outtree->Branch("gg_DR", &t->gg_DR, "gg_DR/F");
    outtree->Branch("mjj_wkinfit", &t->mjj_wkinfit, "mjj_wkinfit/F");
    outtree->Branch("mjj_wokinfit", &t->mjj_wokinfit, "mjj_wokinfit/F");
    outtree->Branch("mtot_wokinfit", &t->mtot_wokinfit, "mtot_wokinfit/F");
    outtree->Branch("cut_based_ct", &t->cut_based_ct, "cut_based_ct/I");
    outtree->Branch("evWeight", &t->evWeight_w_btagSF, "evWeight/F");
    outtree->Branch("weight", &t->weight, "weight/F");
    outtree->Branch("weightBtagSF", &t->weightBtagSF, "weightBtagSF/F");
    outtree->Branch("weightBtagSFerrUp", &t->weightBtagSFerrUp, "weightBtagSFerrUp/F");
    outtree->Branch("weightBtagSFerrDown", &t->weightBtagSFerrDown, "weightBtagSFerrDown/F");
    return;
}
