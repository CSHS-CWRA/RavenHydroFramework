/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2018 the Raven Development Team, Ayman Khedr, Konhee Lee
  ----------------------------------------------------------------*/

#include "TimeSeriesABC.h"
#include "Diagnostics.h"

/*****************************************************************
Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/
//////////////////////////////////////////////////////////////////
/// \brief Implementation of the CDiagnostic constructor
/// \param typ [in] type of diagnostics
//
CDiagnostic::CDiagnostic(diag_type typ)
{
  _type =typ;
  _width =DOESNT_EXIST;
}
//////////////////////////////////////////////////////////////////
/// \brief Implementation of the CDiagnostic constructor
/// \param typ [in] type of diagnostics
//
CDiagnostic::CDiagnostic(diag_type typ, int wid)
{
  _type =typ;
  _width =wid;
}
//////////////////////////////////////////////////////////////////
/// \brief Implementation of the CDiagnostic destructor
//
CDiagnostic::~CDiagnostic(){}
//////////////////////////////////////////////////////////////////
/// \brief returns the name of the diagnostic
//
string CDiagnostic::GetName() const
{
  switch (_type)
  {
  case(DIAG_NASH_SUTCLIFFE):    {return "DIAG_NASH_SUTCLIFFE"; break;}
  case(DIAG_RMSE):              {return "DIAG_RMSE";break;}
  case(DIAG_PCT_BIAS):          {return "DIAG_PCT_BIAS";break;}
  case(DIAG_ABSERR):            {return "DIAG_ABSERR";break;}
  case(DIAG_ABSMAX):            {return "DIAG_ABSMAX";break;}
  case(DIAG_PDIFF):             {return "DIAG_PDIFF";break;}
  case(DIAG_TMVOL):             {return "DIAG_TMVOL";break;}
  case(DIAG_RCOEF):             {return "DIAG_RCOEF"; break;}
  case(DIAG_NSC):               {return "DIAG_NSC";break;}
  case(DIAG_RSR):               {return "DIAG_RSR";break;}
  case(DIAG_R2):                {return "DIAG_R2";break;}
  case(DIAG_CUMUL_FLOW):        {return "DIAG_CUMUL_FLOW";break;}
  case(DIAG_LOG_NASH):          {return "DIAG_LOG_NASH";break;}
  case(DIAG_KLING_GUPTA):       {return "DIAG_KLING_GUPTA";break;}
  case(DIAG_NASH_SUTCLIFFE_DER):{return "DIAG_NASH_SUTCLIFFE_DER"; break;}
  case(DIAG_RMSE_DER):          {return "DIAG_RMSE_DER"; break;}
  case(DIAG_KLING_GUPTA_DER):   {return "DIAG_KLING_GUPTA_DER"; break;}
  case(DIAG_NASH_SUTCLIFFE_RUN):{return"DIAG_NASH_SUTCLIFFE_RUN"; break;}
  default:                      {return "";break;}
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the CDiagnostic constructor
/// \param typ [in] type of diagnostics
//
double CDiagnostic::CalculateDiagnostic(CTimeSeriesABC *pTSMod,
                                        CTimeSeriesABC *pTSObs,
                                        CTimeSeriesABC *pTSWeights,
                                        const optStruct &Options) const
{
  int nn;
  double N=0;
  int width = _width;
  int nVals=pTSObs->GetNumSampledValues();
  double obsval,modval;
  double weight=1;
  int skip=0;
  string filename =pTSObs->GetSourceFile();
  if (!strcmp(pTSObs->GetName().c_str(), "HYDROGRAPH") && (Options.ave_hydrograph == true)){ skip = 1; }//skips first point (initial conditions)
  double dt = Options.timestep;

  switch (_type)
  {
  case(DIAG_NASH_SUTCLIFFE)://-----------------------------------------
  case(DIAG_NASH_SUTCLIFFE_DER):
  case(DIAG_NASH_SUTCLIFFE_RUN):
  {
    if(_type==DIAG_NASH_SUTCLIFFE_DER)
    {
      nVals -= 1;     // Reduce nvals by 1 for derivative of NSE
    }
		else if(_type==DIAG_NASH_SUTCLIFFE_RUN)
		{
			if (width < 2)
			{
				string warn = "Provide average width greater than 1 in format: DIAG_NASH_SUTCLIFFE_RUN[n]";
				WriteWarning(warn, Options.noisy);
				return -ALMOST_INF;
			}
			if (width * 2 > nVals)
			{
				string warn = "Not enough sample values. Check width and timeseries";
				WriteWarning(warn, Options.noisy);
				return -ALMOST_INF;
			}
			nVals -= width;
			skip = width;

		}
   
    bool ValidObs = false;
    bool ValidMod = false;
    bool ValidWeight = false;
    double avg=0;
    N=0;

    for (nn=skip;nn<nVals;nn++)
    {
      if(_type==DIAG_NASH_SUTCLIFFE)
      {
        weight=1.0;
        obsval = pTSObs->GetSampledValue(nn);
        modval = pTSMod->GetSampledValue(nn);
        if(pTSWeights != NULL){ weight=pTSWeights->GetSampledValue(nn);}
        if (obsval != RAV_BLANK_DATA) { ValidObs = true; }
        if (modval != RAV_BLANK_DATA) { ValidMod = true; }
      }
      else if(_type==DIAG_NASH_SUTCLIFFE_DER)
      {
        // obsval and modval becomes (dS(n+1)-dS(n))/dt
        weight=1.0;
        obsval = pTSObs->GetSampledValue(nn+1) - pTSObs->GetSampledValue(nn);
        modval = pTSMod->GetSampledValue(nn+1) - pTSMod->GetSampledValue(nn);
        obsval /= dt;
        modval /= dt;

        // both values at current timestep and next time step need to be valid
        if (pTSObs->GetSampledValue(nn) != RAV_BLANK_DATA && pTSObs->GetSampledValue(nn+1) != RAV_BLANK_DATA) { ValidObs = true; }
        if (pTSMod->GetSampledValue(nn) != RAV_BLANK_DATA && pTSMod->GetSampledValue(nn+1) != RAV_BLANK_DATA) { ValidMod = true; }
        if(pTSWeights != NULL){ weight=(pTSWeights->GetSampledValue(nn+1))*(pTSWeights->GetSampledValue(nn));}
      }
		  else if (_type == DIAG_NASH_SUTCLIFFE_RUN)
		  {
			  int front = 0;
			  int back = 0;
			  double modavg = 0.0;
			  double obsavg = 0.0;
			  weight = 1.0;
        front = (int)(floor(width / 2));
			  if (width % 2 == 1) {
				  back = front;
			  }
			  else{
				  back = front - 1;
			  }

			  for (int k = nn - front; k <= nn + back; k++)
			  {
				  modavg += pTSMod->GetSampledValue(k);
				  obsavg += pTSObs->GetSampledValue(k);
				  if (pTSWeights != NULL) { weight *= pTSWeights->GetSampledValue(k); }
			  }

			  modval = modavg / width;
			  obsval = obsavg / width;
		  }

      if (weight != 0) { ValidWeight = true; }

      if ((obsval != RAV_BLANK_DATA) && (modval != RAV_BLANK_DATA) && weight != 0) {
        avg+=obsval*weight;
        N+= weight;
      }
    }
    avg/=N;

    double sum1(0.0),sum2(0.0);
    for (nn=skip;nn<nVals;nn++)
    {
      if(_type==DIAG_NASH_SUTCLIFFE)
      {
        weight=1.0;
        obsval = pTSObs->GetSampledValue(nn);
        modval = pTSMod->GetSampledValue(nn);
        if(pTSWeights != NULL){ weight=pTSWeights->GetSampledValue(nn);}
        if (weight != 0) { ValidWeight = true; }

        if (obsval != RAV_BLANK_DATA) { ValidObs = true; }
        else {weight = 0;}

        if (modval != RAV_BLANK_DATA) { ValidMod = true; }
        else {weight = 0;}
      }
      else if(_type==DIAG_NASH_SUTCLIFFE_DER)
      {
        // obsval and modval becomes (dS(n+1)-dS(n))/dt
        weight=1.0;
        obsval = pTSObs->GetSampledValue(nn+1) - pTSObs->GetSampledValue(nn);
        modval = pTSMod->GetSampledValue(nn+1) - pTSMod->GetSampledValue(nn);
        obsval /= dt;
        modval /= dt;

        // both values at current timestep and next time step need to be valid
        if(pTSWeights != NULL){ weight=0.5*(pTSWeights->GetSampledValue(nn+1))+(pTSWeights->GetSampledValue(nn));}
        if (weight != 0) { ValidWeight = true; }

        if (pTSObs->GetSampledValue(nn) != RAV_BLANK_DATA && pTSObs->GetSampledValue(nn+1) != RAV_BLANK_DATA) { ValidObs = true;}
        else {weight = 0;}

        if (pTSMod->GetSampledValue(nn) != RAV_BLANK_DATA && pTSMod->GetSampledValue(nn+1) != RAV_BLANK_DATA) { ValidMod = true;}
        else {weight = 0;}
      }
		  else if (_type == DIAG_NASH_SUTCLIFFE_RUN)
		  {
			  int front = 0;
			  int back = 0;
			  double modavg = 0.0;
			  double obsavg = 0.0;
			  weight = 1.0;

        front = (int) floor(width / 2);
			  if (width % 2 == 1) {
				  back = front;
			  }
			  else {
				  back = front - 1;
			  }

			  for (int k = nn - front; k <= nn + back; k++)
			  {
				  modavg += pTSMod->GetSampledValue(k);
				  obsavg += pTSObs->GetSampledValue(k);
				  if (pTSWeights != NULL) { weight *= pTSWeights->GetSampledValue(k); }
			  }

			  modval = modavg / width;
			  obsval = obsavg / width;


			  // Track that all values in width are valid to be used in calculation
			  int validtrack_obs = 0;
			  int validtrack_mod = 0;

			  for (int k = nn - front; k <= nn + back; k++)
			  {
				  if (pTSObs->GetSampledValue(k) != RAV_BLANK_DATA) { validtrack_obs += 1; }
				  if (pTSMod->GetSampledValue(k) != RAV_BLANK_DATA) { validtrack_mod += 1; }
				  if (pTSWeights != NULL) { weight *= pTSWeights->GetSampledValue(k);}
			  }

			  if (validtrack_obs == width) { ValidObs = true; } else{weight = 0;}
			  if (validtrack_mod == width) { ValidMod = true; } else{weight = 0;}
			  if (weight != 0) { ValidWeight = true; }
		
		  }

      sum1 += pow(obsval - modval,2)*weight;
      sum2 += pow(obsval - avg, 2)*weight;
    }


    if (ValidObs && ValidMod && ValidWeight)
    {
      return 1.0 - sum1 / sum2;
    }
    else
    {
      string p1 = "OBSERVED_DATA  "+filename;
      string p2 = "MODELED_DATA  ";
      string p3 = "WEIGHTS ";
      if (ValidObs) { p1 = ""; }
      if (ValidMod) { p2 = ""; }
      if (ValidWeight) { p3 = ""; }

      switch (_type)
      {
      case(DIAG_NASH_SUTCLIFFE):
      {
        string warn = "DIAG_NASH_SUTCLIFFE not performed correctly. Check " + p1 + p2 + p3;
        WriteWarning(warn, Options.noisy);
        break;
      }
      case(DIAG_NASH_SUTCLIFFE_DER):
      {
        string warn = "DIAG_NASH_SUTCLIFFE_DER not performed correctly. Check " + p1 + p2 + p3;
        WriteWarning(warn, Options.noisy);
        break;
      }
		  case(DIAG_NASH_SUTCLIFFE_RUN):
		  {
			  string warn = "DIAG_NASH_SUTCLIFFE_RUN not performed correctly. Check " + p1 + p2 + p3;
			  WriteWarning(warn, Options.noisy);
			  break;
		  }
      }
      return -ALMOST_INF;
    }
    break;
  }
  case(DIAG_RMSE)://----------------------------------------------------
  case(DIAG_RMSE_DER)://----------------------------------------------------
  {
    bool ValidObs = false;
    bool ValidMod = false;
    bool ValidWeight = false;

    double sum;
    N=0;
    sum=0;

    for (nn=skip;nn<nVals;nn++)
    {
      switch (_type)
      {
      case(DIAG_RMSE):
      {
        weight=1.0;
        obsval = pTSObs->GetSampledValue(nn);
        modval = pTSMod->GetSampledValue(nn);
        if(pTSWeights != NULL){ weight=pTSWeights->GetSampledValue(nn);}
        if (weight != 0) { ValidWeight = true; }

        if (obsval != RAV_BLANK_DATA) { ValidObs = true; }
        else {weight = 0;}

        if (modval != RAV_BLANK_DATA) { ValidMod = true; }
        else {weight = 0;}

        break;
      }
      case(DIAG_RMSE_DER):
      {
        weight=1.0;
        obsval = pTSObs->GetSampledValue(nn+1) - pTSObs->GetSampledValue(nn);
        modval = pTSMod->GetSampledValue(nn+1) - pTSMod->GetSampledValue(nn);
        obsval /= dt;
        modval /= dt;


        // both values at current timestep and next time step need to be valid
        if(pTSWeights != NULL){ weight=(pTSWeights->GetSampledValue(nn+1))*(pTSWeights->GetSampledValue(nn));}
        if (weight != 0) { ValidWeight = true; }

        if (pTSObs->GetSampledValue(nn) != RAV_BLANK_DATA && pTSObs->GetSampledValue(nn+1) != RAV_BLANK_DATA) { ValidObs = true;}
        else {weight = 0;}

        if (pTSMod->GetSampledValue(nn) != RAV_BLANK_DATA && pTSMod->GetSampledValue(nn+1) != RAV_BLANK_DATA) { ValidMod = true;}
        else {weight = 0;}

        break;
      }
      }

      sum+=pow(obsval-modval,2)*weight;
      N+=weight;
    }
    if (ValidObs && ValidMod && ValidWeight) {
      return sqrt(sum / N);
    }
    else
    {
      string p1 = "OBSERVED_DATA  ";
      string p2 = "MODELED_DATA  ";
      string p3 = "WEIGHTS ";
      if (ValidObs) { p1 = ""; }
      if (ValidMod) { p2 = ""; }
      if (ValidWeight) { p3 = ""; }

      switch (_type)
      {
      case(DIAG_RMSE):
      {
        string warn = "DIA_RMSE not performed correctly. Check " + p1 + p2 + p3;
        WriteWarning(warn, Options.noisy);
        break;
      }
      case(DIAG_RMSE_DER):
      {
        string warn = "DIAG_RMSE_DER not performed correctly. Check " + p1 + p2 + p3;
        WriteWarning(warn, Options.noisy);
        break;
      }
      }

      return -ALMOST_INF;
    }
    break;
  }
  case(DIAG_PCT_BIAS)://-------------------------------------------------
  {
    bool ValidObs = false;
    bool ValidMod = false;
    bool ValidWeight = false;

    double sum1(0.0),sum2(0.0);
    N=0;

    for (nn=skip;nn<nVals;nn++)
    {
      obsval=pTSObs->GetSampledValue(nn);
      modval=pTSMod->GetSampledValue(nn);
      weight=1.0;
      if(pTSWeights != NULL){ weight=pTSWeights->GetSampledValue(nn);}

      if (obsval != RAV_BLANK_DATA) { ValidObs = true; }
      if (modval != RAV_BLANK_DATA) { ValidMod = true; }
      if (weight != 0) { ValidWeight = true; }

      if (obsval != RAV_BLANK_DATA && modval != RAV_BLANK_DATA && weight != 0) {
        sum1+=(modval-obsval)*weight;
        sum2+=obsval*weight;
      }
    }
    if (ValidObs && ValidMod && ValidWeight)
    {
      return 100.0*sum1/sum2;
    }
    else
    {
      string p1 = "OBSERVED_DATA  "+filename;
      string p2 = "MODELED_DATA  ";
      string p3 = "WEIGHTS ";
      if (ValidObs) { p1 = ""; }
      if (ValidMod) { p2 = ""; }
      if (ValidWeight) { p3 = ""; }

      string warn = "DIAG_PCT_BIAS not performed correctly. Check " + p1 + p2 + p3;
      WriteWarning(warn, Options.noisy);
      return -ALMOST_INF;
    }
    break;
  }
  case(DIAG_ABSERR) ://----------------------------------------------------
  {
    bool ValidObs = false;
    bool ValidMod = false;
    bool ValidWeight = false;

    double sum;
    N = 0;
    sum = 0.0;

    for (nn = skip; nn < nVals; nn++)
    {
      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      weight=1.0;
      if (pTSWeights != NULL) { weight = pTSWeights->GetSampledValue(nn); }

      if (obsval != RAV_BLANK_DATA) { ValidObs = true; }
      if (modval != RAV_BLANK_DATA) { ValidMod = true; }
      if (weight != 0) { ValidWeight = true; }

      if (obsval != RAV_BLANK_DATA && modval != RAV_BLANK_DATA && weight != 0){
        sum += abs(obsval - modval)*weight;
        N+=weight;
      }
    }

    if (ValidObs && ValidMod && ValidWeight)
    {
      return sum / N;
    }
    else
    {
      string p1 = "OBSERVED_DATA  "+filename;
      string p2 = "MODELED_DATA  ";
      string p3 = "WEIGHTS ";
      if (ValidObs) { p1 = ""; }
      if (ValidMod) { p2 = ""; }
      if (ValidWeight) { p3 = ""; }

      string warn = "DIAG_ABSERR not performed correctly. Check " + p1 + p2 + p3;
      WriteWarning(warn, Options.noisy);
      return -ALMOST_INF;
    }
    break;
  }
  case(DIAG_PDIFF) ://----------------------------------------------------
  {
    bool ValidObs = false;
    bool ValidMod = false;
    bool ValidWeight = false;

    double maxObs = 0;
    double maxMod = 0;

    for (nn = skip; nn<nVals; nn++)
    {
      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      weight=1.0;
      if (pTSWeights != NULL) { weight = pTSWeights->GetSampledValue(nn); }

      if (obsval != RAV_BLANK_DATA) { ValidObs = true; }
      if (modval != RAV_BLANK_DATA) { ValidMod = true; }
      if (weight != 0) { ValidWeight = true; }

      if (obsval != RAV_BLANK_DATA && modval != RAV_BLANK_DATA && weight != 0) {
        if (obsval > maxObs) {
          maxObs = obsval;
        }
        if (modval > maxMod) {
          maxMod = modval;
        }
      }
    }
    if (ValidObs && ValidMod && ValidWeight)
    {
      return maxMod - maxObs;
    }
    else
    {
      string p1 = "OBSERVED_DATA  "+filename;
      string p2 = "MODELED_DATA  ";
      string p3 = "WEIGHTS ";
      if (ValidObs) { p1 = ""; }
      if (ValidMod) { p2 = ""; }
      if (ValidWeight) { p3 = ""; }

      string warn = "DIAG_PDIFF not performed correctly. Check " + p1 + p2 + p3;
      WriteWarning(warn, Options.noisy);
      return -ALMOST_INF;
    }
    break;
  }
  case(DIAG_ABSMAX) ://----------------------------------------------------
  {
    bool ValidObs = false;
    bool ValidMod = false;
    bool ValidWeight = false;

    double maxerr = 0;

    for (nn = skip; nn<nVals; nn++)
    {
      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      weight=1.0;
      if (pTSWeights != NULL) { weight = pTSWeights->GetSampledValue(nn); }

      if (obsval != RAV_BLANK_DATA) { ValidObs = true; }
      if (modval != RAV_BLANK_DATA) { ValidMod = true; }
      if (weight != 0) { ValidWeight = true; }

      if (obsval != RAV_BLANK_DATA && modval != RAV_BLANK_DATA && weight != 0) {
        if (abs(obsval - modval) > maxerr ) {
          maxerr = abs(obsval - modval);
        }
      }
    }
    if (ValidObs && ValidMod && ValidWeight)
    {
      return maxerr;
    }
    else
    {
      string p1 = "OBSERVED_DATA  "+filename;
      string p2 = "MODELED_DATA  ";
      string p3 = "WEIGHTS ";
      if (ValidObs) { p1 = ""; }
      if (ValidMod) { p2 = ""; }
      if (ValidWeight) { p3 = ""; }

      string warn = "DIAG_ABSMAX not performed correctly. Check " + p1 + p2 + p3;
      WriteWarning(warn, Options.noisy);
      return -ALMOST_INF;
    }
    break;
  }
  case(DIAG_TMVOL) ://----------------------------------------------------
  {
    bool ValidObs = false;
    bool ValidMod = false;
    bool ValidWeight = false;

    int mon = 0;        // Current Month
    double n_days = 0;     // Number of days at the current month
    double tmvol = 0;   // Total Monthly Mean Error
    double tempsum = 0; // Temporary sum of errors in a month

    // Find month of first valid entry
    for (nn = skip; nn < nVals; nn++)
    {
      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      weight=1.0;
      if (pTSWeights != NULL) { weight = pTSWeights->GetSampledValue(nn); }

      if (obsval != RAV_BLANK_DATA && modval != RAV_BLANK_DATA && weight != 0)
      {
        time_struct tt;
        JulianConvert(nn, Options.julian_start_day, Options.julian_start_year, tt);
        mon = tt.month;
        break;
      }
    }

    // Perform diagnostics
    for (nn = skip; nn < nVals; nn++)
    {
      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      weight=1.0;
      if (pTSWeights != NULL) { weight = pTSWeights->GetSampledValue(nn); }

      if (obsval != RAV_BLANK_DATA) { ValidObs = true; }
      if (modval != RAV_BLANK_DATA) { ValidMod = true; }
      if (weight != 0) { ValidWeight = true; }

      if (obsval != RAV_BLANK_DATA && modval != RAV_BLANK_DATA && weight != 0)
      {
        time_struct tt2;
        JulianConvert(nn, Options.julian_start_day, Options.julian_start_year, tt2);
        // When in the same month
        if (tt2.month == mon)
        {
          tempsum += (modval - obsval)*weight;
          n_days+=weight;
        }
        else
        {
          // Change the current month
          mon = tt2.month;
          // Add up the TMVOL of last month
          if (n_days > 0)
          {
            tmvol += pow((tempsum / n_days), 2);
          }
          // Update tempsum
          tempsum = (modval - obsval)*weight;
          // Update n_days
          n_days = weight;
        }
      }
    }
    // Add up the last month
    if (n_days > 0)
    {
      tmvol += pow(tempsum / n_days, 2);
    }

    if (ValidObs && ValidMod && ValidWeight)
    {
      return tmvol;
    }
    else
    {
      string p1 = "OBSERVED_DATA  "+filename;
      string p2 = "MODELED_DATA  ";
      string p3 = "WEIGHTS ";
      if (ValidObs) { p1 = ""; }
      if (ValidMod) { p2 = ""; }
      if (ValidWeight) { p3 = ""; }

      string warn = "DIAG_TMVOL not performed correctly. Check " + p1 + p2 + p3;
      WriteWarning(warn, Options.noisy);
      return -ALMOST_INF;
    }
    break;
  }
  case(DIAG_RCOEF) ://----------------------------------------------------
  {
    bool ValidObs = false;
    bool ValidMod = false;
    bool ValidWeight = false;

    double ModSum = 0;
    double ObsSum = 0;
    double TopSum = 0;

    for (nn = skip; nn < nVals - 1; nn++)
    {
      double nxtobsval = pTSObs->GetSampledValue(nn + 1);
      double nxtmodval = pTSMod->GetSampledValue(nn + 1);
      double nxtweight=1;
      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      weight=1.0;
      if (pTSWeights != NULL) { weight = pTSWeights->GetSampledValue(nn); nxtweight = pTSWeights->GetSampledValue(nn + 1); }

      if (obsval != RAV_BLANK_DATA && nxtobsval != RAV_BLANK_DATA) { ValidObs = true; }
      if (modval != RAV_BLANK_DATA && nxtmodval != RAV_BLANK_DATA) { ValidMod = true; }
      if (weight != 0 && nxtweight != 0) { ValidWeight = true; }

      if (obsval != RAV_BLANK_DATA && modval != RAV_BLANK_DATA && nxtobsval != RAV_BLANK_DATA && nxtmodval != RAV_BLANK_DATA && weight != 0 && nxtweight != 0)
      {
        ModSum += modval;
        ObsSum += obsval;
        ++N;
      }
    }

    double ModAvg = ModSum / N;
    double ObsAvg = ObsSum / N;

    double ModDiffSum = 0;
    double ObsDiffSum = 0;

    for (nn = skip; nn < nVals - 1; nn++)
    {

      double nxtobsval = pTSObs->GetSampledValue(nn + 1);
      double nxtmodval = pTSMod->GetSampledValue(nn + 1);
      double nxtweight=1.0;
      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      weight=1.0;
      if (pTSWeights != NULL) { weight = pTSWeights->GetSampledValue(nn); nxtweight = pTSWeights->GetSampledValue(nn + 1); }

      if (obsval != RAV_BLANK_DATA && nxtobsval != RAV_BLANK_DATA) { ValidObs = true; }
      if (modval != RAV_BLANK_DATA && nxtmodval != RAV_BLANK_DATA) { ValidMod = true; }
      if (weight != 0 && nxtweight != 0) { ValidWeight = true; }

      if (obsval != RAV_BLANK_DATA && modval != RAV_BLANK_DATA && nxtobsval != RAV_BLANK_DATA && nxtmodval != RAV_BLANK_DATA && weight != 0 && nxtweight != 0)
      {
        TopSum += (nxtmodval - nxtobsval) * (modval - obsval);
        ModDiffSum += pow((modval - ModAvg), 2);
        ObsDiffSum += pow((obsval - ObsAvg), 2);
      }
    }

    double modstd = sqrt((ModDiffSum / N));
    double obsstd = sqrt((ObsDiffSum / N));

    if (ValidObs && ValidMod && ValidWeight)
    {
      return TopSum / N / ((modstd)* (obsstd));
    }
    else
    {
      string p1 = "OBSERVED_DATA  "+filename;
      string p2 = "MODELED_DATA  ";
      string p3 = "WEIGHTS ";
      if (ValidObs) { p1 = ""; }
      if (ValidMod) { p2 = ""; }
      if (ValidWeight) { p3 = ""; }

      string warn = "DIAG_RCOEF not performed correctly. Check " + p1 + p2 + p3;
      WriteWarning(warn, Options.noisy);
      return -ALMOST_INF;
    }
    break;

  }
  case(DIAG_NSC) ://----------------------------------------------------
  {
    bool ValidObs = false;
    bool ValidMod = false;
    bool ValidWeight = false;
    // Counting number of sign changes
    double nsc = 0;

    for (nn = skip; nn < nVals - 1; nn++)
    {

      double nxtobsval = pTSObs->GetSampledValue(nn + 1);
      double nxtmodval = pTSMod->GetSampledValue(nn + 1);
      double nxtweight=1.0;
      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      weight=1.0;
      if (pTSWeights != NULL) { weight = pTSWeights->GetSampledValue(nn); nxtweight = pTSWeights->GetSampledValue(nn + 1); }

      if (obsval != RAV_BLANK_DATA && nxtobsval != RAV_BLANK_DATA) { ValidObs = true; }
      if (modval != RAV_BLANK_DATA && nxtmodval != RAV_BLANK_DATA) { ValidMod = true; }
      if (weight != 0 && nxtweight != 0) { ValidWeight = true; }

      if (obsval != RAV_BLANK_DATA && modval != RAV_BLANK_DATA && nxtobsval != RAV_BLANK_DATA && nxtmodval != RAV_BLANK_DATA && weight != 0 && nxtweight != 0)
      {
        if ((ceil((obsval - modval)*1000)/1000)*(ceil((nxtobsval - nxtmodval)*1000)/1000) < 0)
        {
          ++nsc;
        }
      }
    }

    if (ValidObs && ValidMod && ValidWeight)
    {
      return nsc;
    }
    else
    {
      string p1 = "OBSERVED_DATA  "+filename;
      string p2 = "MODELED_DATA  ";
      string p3 = "WEIGHTS ";
      if (ValidObs) { p1 = ""; }
      if (ValidMod) { p2 = ""; }
      if (ValidWeight) { p3 = ""; }

      string warn = "DIAG_NSC not performed correctly. Check " + p1 + p2 + p3;
      WriteWarning(warn, Options.noisy);
      return -ALMOST_INF;
    }
    return nsc;
    break;

  }
  case(DIAG_RSR) ://----------------------------------------------------
  {
    bool ValidObs = false;
    bool ValidMod = false;
    bool ValidWeight = false;

    double ObsSum = 0;
    double TopSum = 0;
    double BotSum = 0;

    for (nn = skip; nn < nVals; nn++)
    {
      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      weight=1.0;
      if (pTSWeights != NULL) { weight = pTSWeights->GetSampledValue(nn); }

      if (obsval != RAV_BLANK_DATA) { ValidObs = true; }
      if (modval != RAV_BLANK_DATA) { ValidMod = true; }
      if (weight != 0) { ValidWeight = true; }

      if (obsval != RAV_BLANK_DATA && modval != RAV_BLANK_DATA && weight != 0)
      {
        ObsSum += obsval*weight;
        N+=weight;
      }
    }

    double ObsAvg = ObsSum / N;

    for (nn = skip; nn < nVals; nn++)
    {

      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      weight=1.0;
      if (pTSWeights != NULL) { weight = pTSWeights->GetSampledValue(nn); }

      if (obsval != RAV_BLANK_DATA) { ValidObs = true; }
      if (modval != RAV_BLANK_DATA) { ValidMod = true; }
      if (weight != 0) { ValidWeight = true; }

      if (obsval != RAV_BLANK_DATA && modval != RAV_BLANK_DATA && weight != 0)
      {
        TopSum += pow((obsval - modval),2)*weight;
        BotSum += pow((obsval - ObsAvg),2)*weight;
      }
    }
    if (ValidObs && ValidMod && ValidWeight)
    {
      return sqrt(TopSum/BotSum);
    }
    else
    {
      string p1 = "OBSERVED_DATA  "+filename;
      string p2 = "MODELED_DATA  ";
      string p3 = "WEIGHTS ";
      if (ValidObs) { p1 = ""; }
      if (ValidMod) { p2 = ""; }
      if (ValidWeight) { p3 = ""; }

      string warn = "DIAG_RSR not performed correctly. Check " + p1 + p2 + p3;
      WriteWarning(warn, Options.noisy);
      return -ALMOST_INF;
    }
    break;

  }
  case(DIAG_R2) ://----------------------------------------------------
  {
    bool ValidObs = false;
    bool ValidMod = false;
    bool ValidWeight = false;

    double ObsSum = 0;
    double ModSum = 0;

    double CovXY = 0;
    double CovXX = 0;
    double CovYY = 0;

    for (nn = skip; nn < nVals; nn++)
    {
      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      weight=1.0;
      if (pTSWeights != NULL) { weight = pTSWeights->GetSampledValue(nn); }

      if (obsval != RAV_BLANK_DATA) { ValidObs = true; }
      if (modval != RAV_BLANK_DATA) { ValidMod = true; }
      if (weight != 0) { ValidWeight = true; }

      if (obsval != RAV_BLANK_DATA && modval != RAV_BLANK_DATA && weight != 0)
      {
        ObsSum += obsval*weight;
        ModSum += modval*weight;
        N+=weight;
      }
    }

    double ObsAvg = ObsSum / N;
    double ModAvg = ModSum / N;

    for (nn = skip; nn < nVals; nn++)
    {

      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      weight=1.0;
      if (pTSWeights != NULL) { weight = pTSWeights->GetSampledValue(nn); }

      if (obsval != RAV_BLANK_DATA) { ValidObs = true; }
      if (modval != RAV_BLANK_DATA) { ValidMod = true; }
      if (weight != 0) { ValidWeight = true; }

      if (obsval != RAV_BLANK_DATA && modval != RAV_BLANK_DATA && weight != 0)
      {
        CovXY += weight*(modval-ModAvg)*(obsval-ObsAvg);
        CovXX += weight*(modval-ModAvg)*(modval-ModAvg);
        CovYY += weight*(obsval-ObsAvg)*(obsval-ObsAvg);
      }
    }

    CovXY /= N;
    CovXX /= N;
    CovYY /= N;

    if (ValidObs && ValidMod && ValidWeight)
    {
      return pow(CovXY,2)/(CovXX*CovYY);
    }
    else
    {
      string p1 = "OBSERVED_DATA  "+filename;
      string p2 = "MODELED_DATA  ";
      string p3 = "WEIGHTS ";
      if (ValidObs) { p1 = ""; }
      if (ValidMod) { p2 = ""; }
      if (ValidWeight) { p3 = ""; }

      string warn = "DIAG_R2 not performed correctly. Check " + p1 + p2 + p3;
      WriteWarning(warn, Options.noisy);
      return -ALMOST_INF;
    }
    break;

  }

  case(DIAG_LOG_NASH)://-----------------------------------------
  {
    bool ValidObs = false;
    bool ValidMod = false;
    bool ValidWeight = false;
    double avg=0;
    N=0;

    for (nn=skip;nn<nVals;nn++)
    {
      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      weight=1.0;
      if (pTSWeights != NULL) { weight = pTSWeights->GetSampledValue(nn); }

      //log transformation
      if (obsval <= 0.0){obsval=RAV_BLANK_DATA;} //negative values treated as invalid observations/modeled
      else              {obsval=log(obsval);               }
      if (modval <= 0.0){modval=RAV_BLANK_DATA;}
      else              {modval=log(modval);               }

      if (obsval != RAV_BLANK_DATA) { ValidObs = true; }
      if (modval != RAV_BLANK_DATA) { ValidMod = true; }
      if (weight != 0) { ValidWeight = true; }

      if ((obsval != RAV_BLANK_DATA) && (modval != RAV_BLANK_DATA) && (weight != 0))
      {
        avg+=obsval*weight;
        N+= weight;
      }
    }
    avg/=N;

    double sum1(0.0),sum2(0.0);
    for (nn=skip;nn<nVals;nn++)
    {
      obsval=pTSObs->GetSampledValue(nn);
      modval=pTSMod->GetSampledValue(nn);
      weight=1.0;
      if(pTSWeights != NULL){ weight=pTSWeights->GetSampledValue(nn);}

      //log transformation
      if (obsval <= 0.0){obsval=RAV_BLANK_DATA;} //negative values treated as invalid observations/modeled
      else              {obsval=log(obsval);               }
      if (modval <= 0.0){modval=RAV_BLANK_DATA;}
      else              {modval=log(modval);               }

      if ((obsval != RAV_BLANK_DATA) && (modval != RAV_BLANK_DATA) && (weight != 0)) {
        sum1+=pow(obsval-modval,2)*weight;
        sum2+=pow(obsval-avg   ,2)*weight;
      }
    }
    if (ValidObs && ValidMod && ValidWeight)
    {
      return 1.0 - sum1 / sum2;
    }
    else
    {
      string p1 = "OBSERVED_DATA  "+filename;
      string p2 = "MODELED_DATA  ";
      string p3 = "WEIGHTS ";
      if (ValidObs) { p1 = ""; }
      if (ValidMod) { p2 = ""; }
      if (ValidWeight) { p3 = ""; }

      string warn = "DIAG_LOG_NASH not performed correctly. Check " + p1 + p2 + p3;
      WriteWarning(warn, Options.noisy);
      return -ALMOST_INF;
    }
    break;
  }
  case(DIAG_CUMUL_FLOW)://-----------------------------------------
  { // % error in cumulative flows
    bool ValidObs = false;
    bool ValidMod = false;
    bool ValidWeight = false;
    double sumobs=0;
    double summod=0;
    N=0;

    for (nn=skip;nn<nVals;nn++)
    {
      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      weight=1.0;
      if (pTSWeights != NULL) { weight = pTSWeights->GetSampledValue(nn); }

      if (obsval != RAV_BLANK_DATA) { ValidObs = true; }
      if (modval != RAV_BLANK_DATA) { ValidMod = true; }
      if (weight != 0) { ValidWeight = true; }

      if ((obsval != RAV_BLANK_DATA) && (modval != RAV_BLANK_DATA) && (weight != 0))
      {
        sumobs+=obsval*weight;
        summod+=modval*weight;
      }
    }

    if (ValidObs && ValidMod && ValidWeight)
    {

      return (summod-sumobs)/sumobs; //relative error in cumulative flow
    }
    else
    {
      string p1 = "OBSERVED_DATA  "+filename;if (ValidObs   ) { p1 = ""; }
      string p2 = "MODELED_DATA  ";          if (ValidMod   ) { p2 = ""; }
      string p3 = "WEIGHTS ";                if (ValidWeight) { p3 = ""; }

      string warn = "DIAG_CUMUL_FLOW not performed correctly. Check " + p1 + p2 + p3;
      WriteWarning(warn, Options.noisy);
      return ALMOST_INF;
    }
    break;
  }
  case(DIAG_KLING_GUPTA)://-----------------------------------------
  case(DIAG_KLING_GUPTA_DER):
  {
    if(_type==DIAG_KLING_GUPTA_DER){ nVals -= 1; }      // Reduce nvals by 1 for derivative of Kling Gupta
    bool ValidObs = false;
    bool ValidMod = false;
    bool ValidWeight = false;
    double ObsAvg = 0;
    double ModAvg = 0;
    N = 0;
    double ObsSum = 0;
    double ModSum = 0;
    double ObsStd = 0;
    double ModStd = 0;
    double Cov = 0;

    for (nn = skip; nn < nVals; nn++)
    {
      weight = 1.0;
      obsval=modval=0.0;
      if (_type==DIAG_KLING_GUPTA)
      {
        obsval = pTSObs->GetSampledValue(nn);
        modval = pTSMod->GetSampledValue(nn);
        if (pTSWeights != NULL) { weight = pTSWeights->GetSampledValue(nn); }
        if (weight != 0) { ValidWeight = true; }

        if (obsval != RAV_BLANK_DATA) { ValidObs = true; }
        else { weight = 0; }

        if (modval != RAV_BLANK_DATA) { ValidMod = true; }
        else { weight = 0; }
      }
      else if (_type==DIAG_KLING_GUPTA_DER)
      {
        // obsval and modval becomes (dS(n+1)-dS(n))/dt
        obsval = pTSObs->GetSampledValue(nn + 1) - pTSObs->GetSampledValue(nn);
        modval = pTSMod->GetSampledValue(nn + 1) - pTSMod->GetSampledValue(nn);
        obsval /= dt;
        modval /= dt;

        // both values at current timestep and next time step need to be valid
        if (pTSWeights != NULL) { weight = (pTSWeights->GetSampledValue(nn + 1))*(pTSWeights->GetSampledValue(nn)); }
        if (weight != 0) { ValidWeight = true; }

        if (pTSObs->GetSampledValue(nn) != RAV_BLANK_DATA && pTSObs->GetSampledValue(nn + 1) != RAV_BLANK_DATA) { ValidObs = true; }
        else { weight = 0; }

        if (pTSMod->GetSampledValue(nn) != RAV_BLANK_DATA && pTSMod->GetSampledValue(nn + 1) != RAV_BLANK_DATA) { ValidMod = true; }
        else { weight = 0; }
      }
      ObsSum += obsval*weight;
      ModSum += modval*weight;
      N += weight;
    }
    ObsAvg = ObsSum / N;
    ModAvg = ModSum / N;

    for (nn = skip; nn < nVals; nn++)
    {
      obsval=modval=0.0;
      if (_type==DIAG_KLING_GUPTA)
      {
        weight = 1.0;
        obsval = pTSObs->GetSampledValue(nn);
        modval = pTSMod->GetSampledValue(nn);
        if (pTSWeights != NULL) { weight = pTSWeights->GetSampledValue(nn); }
        if (weight != 0) { ValidWeight = true; }

        if (obsval != RAV_BLANK_DATA) { ValidObs = true; }
        else { weight = 0; }

        if (modval != RAV_BLANK_DATA) { ValidMod = true; }
        else { weight = 0; }
      }
      else if (_type==DIAG_KLING_GUPTA_DER)
      {
        // obsval and modval becomes (dS(n+1)-dS(n))/dt
        weight = 1.0;
        obsval = pTSObs->GetSampledValue(nn + 1) - pTSObs->GetSampledValue(nn);
        modval = pTSMod->GetSampledValue(nn + 1) - pTSMod->GetSampledValue(nn);
        obsval /= dt;
        modval /= dt;


        // both values at current timestep and next time step need to be valid
        if (pTSWeights != NULL) { weight = (pTSWeights->GetSampledValue(nn + 1))*(pTSWeights->GetSampledValue(nn)); }
        if (weight != 0) { ValidWeight = true; }

        if (pTSObs->GetSampledValue(nn) != RAV_BLANK_DATA && pTSObs->GetSampledValue(nn + 1) != RAV_BLANK_DATA) { ValidObs = true; }
        else { weight = 0; }

        if (pTSMod->GetSampledValue(nn) != RAV_BLANK_DATA && pTSMod->GetSampledValue(nn + 1) != RAV_BLANK_DATA) { ValidMod = true; }
        else { weight = 0; }
      }
      ObsStd += pow((obsval - ObsAvg), 2)*weight;
      ModStd += pow((modval - ModAvg), 2)*weight;
      Cov += (obsval - ObsAvg) * (modval - ModAvg)*weight;
    }


    ObsStd = sqrt(ObsStd / N);   // Standard Deviation for Observed Flow
    ModStd = sqrt(ModStd / N);   // Standard Deviation for Modelled Flow
    Cov /= N;                    // Covariance between observed and modelled flows

    double r = Cov / ObsStd / ModStd; // pearson product-moment correlation coefficient
    double Beta = ModAvg / ObsAvg;
    double Alpha = ModStd / ObsStd;

    if (ValidObs && ValidMod && ValidWeight)
    {
      return 1 - sqrt(pow((r - 1), 2) + pow((Alpha - 1), 2) + pow((Beta - 1), 2));
    }
    else
    {
      string p1 = "OBSERVED_DATA  "+filename;
      string p2 = "MODELED_DATA  ";
      string p3 = "WEIGHTS ";
      if (ValidObs) { p1 = ""; }
      if (ValidMod) { p2 = ""; }
      if (ValidWeight) { p3 = ""; }

      string warn;
      if (_type==DIAG_KLING_GUPTA){warn = "DIAG_KLING_GUPTA not performed correctly. Check " + p1 + p2 + p3;}
      else if (_type==DIAG_KLING_GUPTA_DER){warn = "DIAG_KLING_GUPTA_DER not performed correctly. Check " + p1 + p2 + p3;}
      WriteWarning(warn, Options.noisy);

      return -ALMOST_INF;
    }
    break;
  }
  default:
  {
    return 0.0; break;
  }
  }//end switch


}
