Delphes_CEPC
=======

Delphes is a C++ framework, performing a fast multipurpose detector response simulation.

More details can be found on the Delphes website http://cp3.irmp.ucl.ac.be/projects/delphes

This is simulation of 4th detector at CEPC with Delphes 

The changes are recorded in CEPC_CHANGELOG

Draw pictures path: drawcode/ 



Quick start with Delphes
========================

   download Delphes_CEPC first
   then
```
   tar -zxf Delphes_CEPC-1.0.tar.gz
```
Configure Delphes
```
source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc10-opt/setup.sh

```
Commands to compile the code:

```
   cd Delphes_CEPC-1.0

   make
```

Finally, we can run Delphes:

```
   ./DelphesSTDHEP
```

Command line parameters:

```
   ./DelphesSTDHEP config_file output_file [input_file(s)]
     config_file - configuration file in Tcl format
     output_file - output file in ROOT format,
     input_file(s) - input file(s) in STDHEP format,
     with no input_file, or when input_file is -, read standard input.
```
some stdhep path:
/cefs/data/stdhep/CEPC240/higgs/update_from_LiangHao_1M/data/


For more detailed documentation, please visit https://cp3.irmp.ucl.ac.be/projects/delphes/wiki/WorkBook


Simple analysis using TTree::Draw
=================================

Now we can start [ROOT](root.cern) and look at the data stored in the output ROOT file.

Start ROOT and load Delphes shared library:

```
   root -l
   gSystem->Load("libDelphes");
```

Open ROOT file and do some basic analysis using Draw or TBrowser:

```
   TFile *f = TFile::Open("delphes_output.root");
   f->Get("Delphes")->Draw("Electron.PT");
   TBrowser browser;
```

Notes:
* ```Delphes``` is the tree name. It can be learned e.g. from TBrowser.
* ```Electron```is the branch name; ```PT``` is a variable (leaf) of this branch.

Complete description of all branches can be found in [doc/RootTreeDescription.html](doc/RootTreeDescription.html).
This information is also available [in the workbook](https://cp3.irmp.ucl.ac.be/projects/delphes/wiki/WorkBook/RootTreeDescription).

Macro-based analysis
====================

Analysis macro consists of histogram booking, event loop (histogram filling),
histogram display.

Start ROOT and load Delphes shared library:

```
   root -l
   gSystem->Load("libDelphes");
```

Basic analysis macro:

```
{
  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add("delphes_output.root");
  
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  // Get pointers to branches used in this analysis
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");

  // Book histograms
  TH1 *histElectronPT = new TH1F("electron pt", "electron P_{T}", 50, 0.0, 100.0);

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
  
    // If event contains at least 1 electron
    if(branchElectron->GetEntries() > 0)
    {
      // Take first electron
      Electron *electron = (Electron*) branchElectron->At(0);
      
      // Plot electron transverse momentum
      histElectronPT->Fill(electron->PT);
      
      // Print electron transverse momentum
      cout << electron->PT << endl;
    }

  }

  // Show resulting histograms
  histElectronPT->Draw();
}
```

More advanced macro-based analysis
==================================

The 'examples' directory contains ROOT macros [Example1.C](examples/Example1.C), [Example2.C](examples/Example2.C) and [Example3.C](examples/Example3.C).

Here are the commands to run these ROOT macros:

```
   root -l
   .X examples/Example1.C("delphes_output.root");
```

or

```
   root -l examples/Example1.C'("delphes_output.root")'
```
