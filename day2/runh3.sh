export OPTIONS="-b --configuration json://config-h3.json"

o2-analysistutorial-o2at-h3-0-vzerotemplateexample ${OPTIONS} | o2-analysis-timestamp ${OPTIONS} | o2-analysis-track-propagation ${OPTIONS} | o2-analysis-event-selection ${OPTIONS} | o2-analysis-lf-lambdakzerobuilder ${OPTIONS} | o2-analysis-pid-tpc ${OPTIONS} | o2-analysis-multiplicity-table ${OPTIONS}
