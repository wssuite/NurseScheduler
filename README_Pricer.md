# Numerical tests with the pricer

New methods have been written for a major revision of the pricer. These methods are all included in the files whose names start with "My" (e.g. MyRCSPP.cpp). The tests are run with the main function that appears in src/PricerMain.cpp and which is called with the binary /bin/pricer.

To call the tests on an instance involving 30 nurses on a 4 weeks horizon, you can call the tests with the following command from the root of the project:
```bin/pricer --dir datasets/ --instance n030w4 --his 0 --weeks 1-2-3-4 --sp-type ROSTER --param paramfiles/default.txt```

Several variants of the pricer can be tested by modifying the following four options in the parameter file (paramfiles/default.txt by default):
sortLabelsOption=0 (or 1)
minimumCostFromSinksOption=0 (or 1)
worstCaseCostOption=1 (or 0)
enumeratedSubPathOption=0 (or 1)

