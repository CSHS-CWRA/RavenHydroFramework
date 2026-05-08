/*----------------------------------------------------------------
  RavenCaller Source Code
  Copyright (c) 2026 the Raven Development Team

  calls Raven.dll - used for testing BMI functionality
  ----------------------------------------------------------------*/
// JRC: look into https://yal.cc/cpp-a-very-tiny-dll/
//
#ifdef _WIN32
#include <Windows.h>
#endif
#include <iostream>
#include <fstream>
#include "D:\James\Software\Raven\src\BMI.h"
using namespace std;

typedef void* (*model_create)();
typedef void (*model_destroy)(void*);

int main()
{
  cout << "Calling Raven from RavenCaller..."<<endl;

#ifdef _WIN32
  const WCHAR* dllFilename=L"D:\\James\\Software\\Raven\\src\\ReleaseDLL\\Raven.dll"; //for debugging
  //const WCHAR* dllFilename=L"Raven.dll"; //if using stationary .dll

  HINSTANCE hRavenBMI = LoadLibrary(dllFilename);
  if (!hRavenBMI) {
    cerr<<"RAVENCALLER ERROR: Unable to load dll."<<endl;
    return 1;
  }

  model_create  bmi_model_create  = (model_create )GetProcAddress(hRavenBMI, "bmi_model_create");
  model_destroy bmi_model_destroy = (model_destroy)GetProcAddress(hRavenBMI, "bmi_model_destroy");
    
  if(bmi_model_create==NULL) {
    cerr<<"RAVENCALLER ERROR: NULL pointer to create function"<<endl;
    FreeLibrary(hRavenBMI);
    return 1;
  }
  
  //Create instance of Raven model within BMI wrapper
  bmixx::Bmi* pRavenInstance=(bmixx::Bmi*)bmi_model_create();
  if (pRavenInstance!=NULL)
  {   
    try {
      //Initialize model
      cout<<"Initializing..."<<endl;
      pRavenInstance->Initialize("config.txt");
    
      //Query model configuration:
      vector<string> inputs =pRavenInstance->GetInputVarNames();
      cout<<"ACCESSIBLE INPUTS: "<<endl;
      for (int i=0; i<inputs.size(); i++){cout<<"-"<<inputs[i]<<endl;}

      vector<string> outputs =pRavenInstance->GetOutputVarNames();
      cout<<"ACCESSIBLE OUTPUTS: "<<endl;
      for (int i=0; i<outputs.size(); i++){cout<<"-"<<outputs[i]<<endl;}

      //Run model - this is where coupling occurs
      double duration=20;  //days   //This duration must be less or equal to that indicated in the .rvi file
      double tstep   =1.0; //days   //This time step is assumed to be the same as that in the .rvi file

      int SoilGrid=pRavenInstance->GetVarGrid("SOIL[2]");
      int FlowGrid=pRavenInstance->GetVarGrid("STREAMFLOW");
      int nHRUs=pRavenInstance->GetGridSize(SoilGrid);
      int nSBs=pRavenInstance->GetGridSize(FlowGrid);

      double *SB_values =new double [nSBs];
      double *HRU_values=new double [nHRUs];

      for (double t=0; t<duration; t+=tstep)
      {
        //override states at start of time step:
        for (int k=0;k<nHRUs;k++){
          HRU_values[k]=20; //mm
        }
        pRavenInstance->SetValue("SOIL[2]",HRU_values);

        //Run Raven for one time step
        pRavenInstance->Update();

        //read outputs at end of time step: 
        pRavenInstance->GetValue("STREAMFLOW",SB_values);
        cout<<"streamflow: "<<SB_values[0]<<endl;
      }
      delete [] SB_values;
      delete [] HRU_values;

      //Finalize model
      pRavenInstance->Finalize();
    } 

    //Exception handling
    catch(const std::exception& e) {
      cout<<" RAVEN BMI DLL ERROR: "<<e.what()<<endl;
      try {
        pRavenInstance->Finalize();
      }
      catch(const std::exception& e) {
       cout<<" RAVEN BMI DLL ERROR (UNABLE TO CLEAN UP MEMORY): "<<e.what()<<endl;
      }
      return 0;
    }
  } 
  else{
    cout<<" RAVENCALLER ERROR: Invalid Raven Instantiation :("<<endl;
  }

  //Clean up
  bmi_model_destroy(pRavenInstance);
  FreeLibrary(hRavenBMI);
#endif
  return 0;
}

