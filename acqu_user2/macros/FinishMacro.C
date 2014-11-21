void FinishMacro()  {

  char command[128];
  TA2GammaDeuterium* treeFile = (TA2GammaDeuterium*)(gAN->GetPhysics());
  treeFile->Print();
  treeFile->ProcessEnd();
  printf("End-of-Run macro executing\n");
  sprintf(command,"kill -9 %d",gSystem->GetPid());
  gSystem->Exec(command);

/*
    TString name;
    printf("\nEnd-of-Run macro executing:\n");

    Char_t* file = "ARHistograms.root";
    TFile f1(file,"RECREATE");
    gROOT->GetList()->Write();
    gROOT->GetList()->Print();
    f1.Close();
    printf("done.\n",file);
    printf("All histograms saved to %s\n",file); 
    gSystem->Exit(0);
*/

}
