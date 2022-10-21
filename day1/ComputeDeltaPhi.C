Double_t ComputeDeltaPhi( Double_t phi1, Double_t phi2) {
  //To be completely sure, use inner products
  Double_t x1, y1, x2, y2;
  x1 = TMath::Cos( phi1 );
  y1 = TMath::Sin( phi1 );
  x2 = TMath::Cos( phi2 );
  y2 = TMath::Sin( phi2 );
  Double_t lInnerProd  = x1*x2 + y1*y2;
  Double_t lVectorProd = x1*y2 - x2*y1;
  
  Double_t lReturnVal = 0;
  if( lVectorProd > 1e-8 ) {
    lReturnVal = TMath::ACos(lInnerProd);
  }
  if( lVectorProd < -1e-8 ) {
    lReturnVal = -TMath::ACos(lInnerProd);
  }
  
  if( lReturnVal < -TMath::Pi()/2. ) {
    lReturnVal += 2.*TMath::Pi();
  }
  
  return lReturnVal;
}

