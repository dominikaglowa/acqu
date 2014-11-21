void PeriodMacro()
{
  // Stuff to do at the end of an analysis run
  // Here all spectra are saved to disk
  printf("Events: %d      Events Accepted: %d\n\n",
	 gAN->GetNEvent(), gAN->GetNEvAnalysed() );
}
