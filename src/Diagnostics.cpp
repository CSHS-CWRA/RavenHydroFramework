/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team, Ayman Khedr, Konhee Lee
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
  case(DIAG_NASH_SUTCLIFFE):    {return "DIAG_NASH_SUTCLIFFE"; }
  case(DIAG_DAILY_NSE):         {return "DIAG_DAILY_NSE"; }
  case(DIAG_RMSE):              {return "DIAG_RMSE";}
  case(DIAG_PCT_BIAS):          {return "DIAG_PCT_BIAS";}
  case(DIAG_ABSERR):            {return "DIAG_ABSERR";}
  case(DIAG_ABSMAX):            {return "DIAG_ABSMAX";}
  case(DIAG_PDIFF):             {return "DIAG_PDIFF";}
  case(DIAG_TMVOL):             {return "DIAG_TMVOL";}
  case(DIAG_RCOEF):             {return "DIAG_RCOEF"; }
  case(DIAG_NSC):               {return "DIAG_NSC";}
  case(DIAG_RSR):               {return "DIAG_RSR";}
  case(DIAG_R2):                {return "DIAG_R2";}
  case(DIAG_CUMUL_FLOW):        {return "DIAG_CUMUL_FLOW";}
  case(DIAG_LOG_NASH):          {return "DIAG_LOG_NASH";}
  case(DIAG_KLING_GUPTA):       {return "DIAG_KLING_GUPTA";}
  case(DIAG_NASH_SUTCLIFFE_DER):{return "DIAG_NASH_SUTCLIFFE_DER"; }
  case(DIAG_RMSE_DER):          {return "DIAG_RMSE_DER"; }
  case(DIAG_KLING_GUPTA_DER):   {return "DIAG_KLING_GUPTA_DER"; }
  case(DIAG_NASH_SUTCLIFFE_RUN):{return "DIAG_NASH_SUTCLIFFE_RUN"; }
  case(DIAG_MBF):               {return "DIAG_MBF"; }
  case(DIAG_R4MS4E):            {return "DIAG_R4MS4E"; }
  case(DIAG_RTRMSE):            {return "DIAG_RTRMSE"; }
  case(DIAG_RABSERR):           {return "DIAG_RABSERR"; }
  case(DIAG_PERSINDEX):         {return "DIAG_PERSINDEX"; }
  case(DIAG_NSE4):              {return "DIAG_NSE4"; }
  default:                      {return "";}
  }
}
//////////////////////////////////////////////////////////////////
/// \brief returns the type of the diagnostic
//
diag_type CDiagnostic::GetType() const
{
  return _type;
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the CDiagnostic constructor
/// \param typ [in] type of diagnostics
//
double CDiagnostic::CalculateDiagnostic(CTimeSeriesABC  *pTSMod,
                                        CTimeSeriesABC  *pTSObs,
                                        CTimeSeriesABC  *pTSWeights,
                                        const double    &starttime,
                                        const double    &endtime,
                                        const optStruct &Options) const
{
  int nn;
  double N=0;
  string filename =pTSObs->GetSourceFile();
  double obsval,modval;
  double weight=1;

  int    skip  =0;
  if (!strcmp(pTSObs->GetName().c_str(), "HYDROGRAPH") && (Options.ave_hydrograph == true)){ skip = 1; }//skips first point (initial conditions)
  double dt = Options.timestep;

  int nnstart=pTSObs->GetTimeIndexFromModelTime(starttime+1.0)+skip; //works for avg. hydrographs - unclear why 1.0 is neccessary
  int nnend  =pTSObs->GetTimeIndexFromModelTime(endtime  +1.0)+1; //+1 is just because below loops expressed w.r.t N, not N-1

  switch (_type)
  {
  case(DIAG_NASH_SUTCLIFFE)://----------------------------------------------------
  {
    double avgobs=0.0;
    N=0;

    for(nn=nnstart;nn<nnend;nn++)
    {
      weight=1.0;
      obsval = pTSObs->GetSampledValue(nn);
      if(pTSWeights != NULL) { weight=pTSWeights->GetSampledValue(nn); }
      if(obsval==RAV_BLANK_DATA) { weight=0.0; }
      avgobs+=weight*obsval;
      N     +=weight;
    }
    if(N>0.0) { avgobs/=N; }

    double sum1(0.0),sum2(0.0);
    for(nn=nnstart;nn<nnend;nn++)
    {
      weight=1.0;
      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      if(pTSWeights != NULL) { weight=pTSWeights->GetSampledValue(nn); }
      if((obsval==RAV_BLANK_DATA) || (modval==RAV_BLANK_DATA)) { weight=0.0; }

      sum1 += weight*pow(obsval - modval,2);
      sum2 += weight*pow(obsval - avgobs,2);
    }

    if(N>0)
    {
      return 1.0 - (sum1 / sum2);
    }
    else
    {
      string warn = "DIAG_NASH_SUTCLIFFE not calculated. Missing non-zero weighted observations during simulation duration.";
      WriteWarning(warn,Options.noisy);
      return -ALMOST_INF;
    }
    break;
  }
  case(DIAG_NASH_SUTCLIFFE_DER)://----------------------------------------------------
  {
    nnend -= 1;     // Reduce nnend by 1 for derivative of NSE
    double obsval2,modval2;
    double avg=0;
    N=0.0;
    
    for(nn=nnstart;nn<nnend;nn++)
    {
      weight=1.0;
      obsval = pTSObs->GetSampledValue(nn);
      obsval2= pTSObs->GetSampledValue(nn+1);
      if(pTSWeights != NULL) { weight=(pTSWeights->GetSampledValue(nn+1))*(pTSWeights->GetSampledValue(nn)); }
      if((obsval == RAV_BLANK_DATA) || (obsval2 == RAV_BLANK_DATA)) { weight=0.0; }
     
      avg+= weight*(obsval2-obsval)/dt;
      N  += weight;
    }
    if(N>0.0) {
      avg/=N;
    }
    double sum1(0.0),sum2(0.0);
    for(nn=nnstart;nn<nnend;nn++)
    {
      weight=1.0;
      obsval2 = pTSObs->GetSampledValue(nn+1);
      obsval  = pTSObs->GetSampledValue(nn);
      modval2 = pTSMod->GetSampledValue(nn+1);
      modval  = pTSMod->GetSampledValue(nn);
      if(pTSWeights != NULL) { weight=0.5*(pTSWeights->GetSampledValue(nn+1))+(pTSWeights->GetSampledValue(nn)); }

      if ((obsval==RAV_BLANK_DATA) || (obsval2==RAV_BLANK_DATA) || (modval==RAV_BLANK_DATA) || (modval2==RAV_BLANK_DATA)){weight=0.0;}
        
      sum1 += weight*pow((obsval2-obsval)/dt - (modval2-modval)/dt,2);
      sum2 += weight*pow((obsval2-obsval)/dt - avg,2);
    }

    if(N>0.0)
    {
      return 1.0 - sum1 / sum2;
    }
    else
    {
      string warn = "DIAG_NASH_SUTCLIFFE_DER not calculated. Missing non-zero weighted observations during simulation duration.";
      WriteWarning(warn,Options.noisy);
      return -ALMOST_INF;
    }
    break;
  }
  case(DIAG_NASH_SUTCLIFFE_RUN)://----------------------------------------------------
  {
    if (_width < 2)
		{
			string warn = "Provide average _width greater than 1 in format: DIAG_NASH_SUTCLIFFE_RUN[n]";
			WriteWarning(warn, Options.noisy);
			return -ALMOST_INF;
		}
		if (_width * 2 > nnend)
		{
			string warn = "Not enough sample values. Check width and timeseries";
			WriteWarning(warn, Options.noisy);
			return -ALMOST_INF;
		}
		nnend -= _width;
		nnstart  += _width;
		
    bool ValidObs = false;
    bool ValidMod = false;
    bool ValidWeight = false;
    double avg=0;
    N=0;

    for (nn=nnstart;nn<nnend;nn++)
    {
    	int front = 0;
			int back = 0;
			double modavg = 0.0;
			double obsavg = 0.0;
			weight = 1.0;
      front = (int)(floor(_width / 2));
			if (_width % 2 == 1) { back = front;  }
			else                 { back = front-1;}

			for (int k = nn - front; k <= nn + back; k++)
			{
				modavg += pTSMod->GetSampledValue(k);
				obsavg += pTSObs->GetSampledValue(k);
				if (pTSWeights != NULL) { weight *= pTSWeights->GetSampledValue(k); }
			}

			modval = modavg / _width;
			obsval = obsavg / _width;
		
      if (weight != 0) { ValidWeight = true; }

      if ((obsval != RAV_BLANK_DATA) && (modval != RAV_BLANK_DATA) && weight != 0) {
        avg+=obsval*weight;
        N  += weight;
      }
    }
    avg/=N;

    double sum1(0.0),sum2(0.0);
    for (nn=nnstart;nn<nnend;nn++)
    {
      int front = 0;
      int back = 0;
      double modavg = 0.0;
      double obsavg = 0.0;
      weight = 1.0;

      front = (int)floor(_width / 2);
      if(_width % 2 == 1) {back = front;}
      else                {back = front - 1;}

      for(int k = nn - front; k <= nn + back; k++)
      {
        modavg += pTSMod->GetSampledValue(k);
        obsavg += pTSObs->GetSampledValue(k);
        if(pTSWeights != NULL) { weight *= pTSWeights->GetSampledValue(k); }
      }

      modval = modavg / _width;
      obsval = obsavg / _width;

      // Track that all values in width are valid to be used in calculation
      int validtrack_obs = 0;
      int validtrack_mod = 0;

      for(int k = nn - front; k <= nn + back; k++)
      {
        if(pTSObs->GetSampledValue(k) != RAV_BLANK_DATA) { validtrack_obs += 1; }
        if(pTSMod->GetSampledValue(k) != RAV_BLANK_DATA) { validtrack_mod += 1; }
        if(pTSWeights != NULL) { weight *= pTSWeights->GetSampledValue(k); }
      }

      if(validtrack_obs == _width) { ValidObs = true; }
      else { weight = 0; }
      if(validtrack_mod == _width) { ValidMod = true; }
      else { weight = 0; }
      if(weight != 0) { ValidWeight = true; }
    
      sum1 += pow(obsval - modval,2)*weight;
      sum2 += pow(obsval - avg   ,2)*weight;
    }

    if(N>0.0)
    {
      return 1.0 - sum1 / sum2;
    }
    else
    {
      string warn = "DIAG_NASH_SUTCLIFFE_RUN not performed correctly. Missing non-zero weighted observations during simulation duration.";
      WriteWarning(warn,Options.noisy);
      return -ALMOST_INF;
    }
    break;
  }
  case(DIAG_DAILY_NSE)://----------------------------------------------------
  {
    double avgobs=0.0;
    int freq=(int)(round(1.0/Options.timestep)); 
    
    int shift=(int)(1.0- Options.julian_start_day-floor(Options.julian_start_day+TIME_CORRECTION))*freq; //no. of timesteps nnstartped at start of simulation (starts on first full day)
    if (nnstart>skip){shift=0;}//using diagnostic period always midnight to midnight, no shift required
    double moddaily=0.0;
    double obsdaily=0.0;
    double dailyN  =0.0;
    double N=0;

    for(nn=nnstart+shift;nn<nnend;nn++) //calculate mean observation value
    {
      weight=1.0;
      obsval = pTSObs->GetSampledValue(nn);
      if(pTSWeights != NULL)    { weight=pTSWeights->GetSampledValue(nn); }
      if(obsval==RAV_BLANK_DATA){ weight=0.0; }
      avgobs+=weight*obsval;
      N     +=weight;
    }
    if(N>0.0) { avgobs/=N; }

    double sum1(0.0),sum2(0.0);
    for(nn=nnstart+shift;nn<nnend;nn++) //calculate numerator and denominator of NSE term
    {
      weight=1.0;
      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      if(pTSWeights != NULL) { weight=pTSWeights->GetSampledValue(nn); }
      if((obsval==RAV_BLANK_DATA) || (modval==RAV_BLANK_DATA)) { weight=0.0; }
      
      obsdaily+=weight*obsval;
      moddaily+=weight*modval;
      dailyN  +=weight;
      if(((nn-nnstart+shift)%freq)==(freq-1)) { //last timestep of day
        if(dailyN>0) {
          obsdaily/=dailyN;
          moddaily/=dailyN;
          sum1 += pow(obsdaily - moddaily,2)*weight;
          sum2 += pow(obsdaily - avgobs  ,2)*weight;
        }
        dailyN=obsdaily=moddaily=0.0; //reset for next day
      }
    }

    if(N>0)
    {
      return 1.0 - (sum1 / sum2);
    }
    else
    {
      string warn = "DIAG_DAILY_NSE not calculated. Missing non-zero weighted observations during simulation duration.";
      WriteWarning(warn,Options.noisy);
      return -ALMOST_INF;
    }
    break;
  }
  case(DIAG_RMSE)://----------------------------------------------------
  {
    double sum=0.0;
    N=0;
    for(nn=nnstart;nn<nnend;nn++)
    {
      weight=1.0;
      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      if(pTSWeights != NULL) { weight=pTSWeights->GetSampledValue(nn); }
      if ((obsval==RAV_BLANK_DATA) || (modval==RAV_BLANK_DATA)){weight=0.0;}
  
      sum+=weight*pow(obsval-modval,2);
      N  +=weight;
    }
    if(N>0.0) {
      return sqrt(sum / N);
    }
    else
    {
      string warn = "DIA_RMSE not not calculated. Missing non-zero weighted observations during simulation duration.";
      WriteWarning(warn,Options.noisy);
      return -ALMOST_INF;
    }
    break;
  }
  case(DIAG_RMSE_DER)://----------------------------------------------------
  {
    double sum;
    N=0;
    sum=0;

    for (nn=nnstart;nn<nnend;nn++)
    {

      weight=1.0;
      obsval = pTSObs->GetSampledValue(nn+1) - pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn+1) - pTSMod->GetSampledValue(nn);
      obsval /= dt;
      modval /= dt;

      // both values at current timestep and next time step need to be valid
      if(pTSWeights != NULL){ weight=(pTSWeights->GetSampledValue(nn+1))*(pTSWeights->GetSampledValue(nn));}

      if (pTSObs->GetSampledValue(nn) != RAV_BLANK_DATA && pTSObs->GetSampledValue(nn+1) != RAV_BLANK_DATA) { }
      else {weight = 0;}

      if (pTSMod->GetSampledValue(nn) != RAV_BLANK_DATA && pTSMod->GetSampledValue(nn+1) != RAV_BLANK_DATA) { }
      else {weight = 0;}

      sum+=weight*pow(obsval-modval,2);
      N  +=weight;
    }
    if (N>0.0) {
      return sqrt(sum / N);
    }
    else
    {
      string warn = "DIAG_RMSE_DER not calculated. Missing non-zero weighted observations during simulation duration.";
      WriteWarning(warn, Options.noisy);
      return -ALMOST_INF;
    }
    break;
  }
  case(DIAG_PCT_BIAS)://-------------------------------------------------
  {
    double sum1(0.0),sum2(0.0);
    N=0;
    for (nn=nnstart;nn<nnend;nn++)
    {
      obsval=pTSObs->GetSampledValue(nn);
      modval=pTSMod->GetSampledValue(nn);
      weight=1.0;
      if(pTSWeights != NULL){ weight=pTSWeights->GetSampledValue(nn);}
      if((obsval==RAV_BLANK_DATA) || (modval==RAV_BLANK_DATA)) { weight=0.0; }

      sum1+=weight*(modval-obsval);
      sum2+=weight*obsval;
      N   +=weight;
    }
    if (N>0.0)
    {
      return 100.0*sum1/sum2;
    }
    else
    {
      string warn = "DIAG_PCT_BIAS not calculated. Missing non-zero weighted observations during simulation duration.";
      WriteWarning(warn, Options.noisy);
      return ALMOST_INF;
    }
    break;
  }
  case(DIAG_ABSERR) ://----------------------------------------------------
  {
    double sum=0.0;
    N = 0;
   
    for (nn = nnstart; nn < nnend; nn++)
    {
      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      weight=1.0;
      if (pTSWeights != NULL) { weight = pTSWeights->GetSampledValue(nn); }
      if((obsval==RAV_BLANK_DATA) || (modval==RAV_BLANK_DATA)) { weight=0.0; }

      sum += weight*fabs(obsval - modval);
      N   += weight;
    }
    if (N>0.0)
    {
      return sum / N;
    }
    else
    {
      string warn = "DIAG_ABSERR not calculated. Missing non-zero weighted observations during simulation duration.";
      WriteWarning(warn, Options.noisy);
      return -ALMOST_INF;
    }
    break;
  }
  case(DIAG_ABSMAX)://----------------------------------------------------
  {
    double maxerr = -ALMOST_INF;
    double N=0.0;

    for(nn = nnstart; nn<nnend; nn++)
    {
      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      weight=1.0;
      if(pTSWeights != NULL) { weight = pTSWeights->GetSampledValue(nn); }
      if((obsval==RAV_BLANK_DATA) || (modval==RAV_BLANK_DATA)) { weight=0.0; }

      if(weight > 0.0) {
        if(fabs(obsval - modval) > maxerr) {
          maxerr = fabs(obsval - modval);
        }
        N+=weight;
      }
    }
    if(N>0.0)
    {
      return maxerr;
    }
    else
    {
      string warn = "DIAG_ABSMAX not calculated. Missing non-zero weighted observations during simulation duration.";
      WriteWarning(warn,Options.noisy);
      return -ALMOST_INF;
    }
    break;
  }
  case(DIAG_PDIFF) ://----------------------------------------------------
  {
    double maxObs = 0;
    double maxMod = 0;
    double N=0;

    for (nn = nnstart; nn<nnend; nn++)
    {
      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      weight=1.0;
      if (pTSWeights != NULL) { weight = pTSWeights->GetSampledValue(nn); }
      if((obsval==RAV_BLANK_DATA) || (modval==RAV_BLANK_DATA)) { weight=0.0; }

      if (weight> 0) {
        if (obsval > maxObs) {maxObs = obsval;}
        if (modval > maxMod) {maxMod = modval;}
        N+=weight;
      }
    }
    if (N>0.0)
    {
      return maxMod - maxObs;
    }
    else
    {
      string warn = "DIAG_PDIFF not calculated. Missing non-zero weighted observations during simulation duration.";
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
    for (nn = nnstart; nn < nnend; nn++)
    {
      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      weight=1.0;
      if (pTSWeights != NULL) { weight = pTSWeights->GetSampledValue(nn); }

      if (obsval != RAV_BLANK_DATA && modval != RAV_BLANK_DATA && weight != 0)
      {
        time_struct tt;
		    JulianConvert(nn, Options.julian_start_day, Options.julian_start_year, Options.calendar, tt);
        mon = tt.month;
        break;
      }
    }

    // Perform diagnostics
    for (nn = nnstart; nn < nnend; nn++)
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
		    JulianConvert(nn, Options.julian_start_day, Options.julian_start_year, Options.calendar, tt2);
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

    for (nn = nnstart; nn < nnend - 1; nn++)
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

    for (nn = nnstart; nn < nnend - 1; nn++)
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
    // Counting number of sign changes of difference between time series
    double nxtobsval,nxtmodval;
    double nsc = 0;
    double N=0;
    for (nn = nnstart; nn < nnend - 1; nn++)
    {

      nxtobsval = pTSObs->GetSampledValue(nn + 1);
      nxtmodval = pTSMod->GetSampledValue(nn + 1);
      obsval    = pTSObs->GetSampledValue(nn);
      modval    = pTSMod->GetSampledValue(nn);
      double nxtweight=1.0;
      weight=1.0;
      if (pTSWeights != NULL) { weight = pTSWeights->GetSampledValue(nn); nxtweight = pTSWeights->GetSampledValue(nn + 1); }
      if((obsval==RAV_BLANK_DATA) || (modval==RAV_BLANK_DATA) || (nxtobsval==RAV_BLANK_DATA) || (nxtmodval==RAV_BLANK_DATA) || (nxtweight==0.0)) { weight=0.0; }

      if (weight>0.0)
      {
        if ((ceil((obsval - modval)*1000)/1000)*(ceil((nxtobsval - nxtmodval)*1000)/1000) < 0)
        {
          nsc++;
        }
      }
      N+=weight;
    }

    if (N>0.0)
    {
      return nsc;
    }
    else
    {
      string warn = "DIAG_NSC  not calculated. Missing non-zero weighted observations during simulation duration.";
      WriteWarning(warn, Options.noisy);
      return ALMOST_INF;
    }
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

    for (nn = nnstart; nn < nnend; nn++)
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

    for (nn = nnstart; nn < nnend; nn++)
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

    for (nn = nnstart; nn < nnend; nn++)
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

    for (nn = nnstart; nn < nnend; nn++)
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

    for (nn=nnstart;nn<nnend;nn++)
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
    for (nn=nnstart;nn<nnend;nn++)
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

    for (nn=nnstart;nn<nnend;nn++)
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
    if(_type==DIAG_KLING_GUPTA_DER){ nnend -= 1; }      // Reduce nnend by 1 for derivative of Kling Gupta
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

    for (nn = nnstart; nn < nnend; nn++)
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

    for (nn = nnstart; nn < nnend; nn++)
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
  case(DIAG_MBF)://----------------------------------------------------
  {
     bool ValidObs = false;
     bool ValidMod = false;
     bool ValidWeight = false;

     double sum;
     N = 0;
     sum = 0.0;

     for(nn = nnstart; nn < nnend; nn++)
     {
       obsval = pTSObs->GetSampledValue(nn);
       modval = pTSMod->GetSampledValue(nn);
       weight=1.0;
       if(pTSWeights != NULL) { weight = pTSWeights->GetSampledValue(nn); }

       if(obsval != RAV_BLANK_DATA) { ValidObs = true; }
       if(modval != RAV_BLANK_DATA) { ValidMod = true; }
       if(weight != 0) { ValidWeight = true; }

       if((obsval != RAV_BLANK_DATA) && (modval != RAV_BLANK_DATA) && (weight != 0)) {
         sum += weight / (1 + pow(((modval-obsval) / static_cast<double>(2.0*obsval)),2));
         N+=weight;
       }
     }

     if(ValidObs && ValidMod && ValidWeight)
     {
       return sum;
     }
     else
     {
       string p1 = "OBSERVED_DATA  "+filename;
       string p2 = "MODELED_DATA  ";
       string p3 = "WEIGHTS ";
       if(ValidObs) { p1 = ""; }
       if(ValidMod) { p2 = ""; }
       if(ValidWeight) { p3 = ""; }

       string warn = "DIAG_ABSERR not performed correctly. Check " + p1 + p2 + p3;
       WriteWarning(warn,Options.noisy);
       return -ALMOST_INF;
     }
     break;
  }
  case(DIAG_R4MS4E)://----------------------------------------------------
  {
    double sum=0.0;
    N=0;
    for(nn=nnstart;nn<nnend;nn++)
    {
      weight=1.0;
      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      if(pTSWeights != NULL) { weight=pTSWeights->GetSampledValue(nn); }
      if ((obsval==RAV_BLANK_DATA) || (modval==RAV_BLANK_DATA)){weight=0.0;}
  
      sum+=weight*pow(obsval-modval,4);
      N  +=weight;
    }
    if(N>0.0) {
      return pow( (sum/N), 0.25); 
    }
    else
    {
      string warn = "DIAG_R4MS4E not not calculated. Missing non-zero weighted observations during simulation duration.";
      WriteWarning(warn,Options.noisy);
      return -ALMOST_INF;
    }
    break;
  }
  case(DIAG_RTRMSE)://----------------------------------------------------
  {
    double sum=0.0;
    N=0;
    for(nn=nnstart;nn<nnend;nn++)
    {
      weight=1.0;
      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      if(pTSWeights != NULL) { weight=pTSWeights->GetSampledValue(nn); }
      if ((obsval==RAV_BLANK_DATA) || (modval==RAV_BLANK_DATA)){weight=0.0;}
  
      sum+=weight*pow( sqrt(obsval)-sqrt(modval), 2);
      N  +=weight;
    }
    if(N>0.0) {
      return sqrt(sum / N);
    }
    else
    {
      string warn = "DIAG_RTRMSE not not calculated. Missing non-zero weighted observations during simulation duration.";
      WriteWarning(warn,Options.noisy);
      return -ALMOST_INF;
    }
    break;
  }
  case(DIAG_RABSERR) ://----------------------------------------------------
  {
    double avgobs=0.0;
    N=0;

    for(nn=nnstart;nn<nnend;nn++)
    {
      weight=1.0;
      obsval = pTSObs->GetSampledValue(nn);
      if(pTSWeights != NULL) { weight=pTSWeights->GetSampledValue(nn); }
      if(obsval==RAV_BLANK_DATA) { weight=0.0; }
      avgobs+=weight*obsval;
      N     +=weight;
    }
    if(N>0.0) { avgobs/=N; }      
      
    double sum1(0.0),sum2(0.0);
    N = 0;
   
    for (nn = nnstart; nn < nnend; nn++)
    {
      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      weight=1.0;
      if (pTSWeights != NULL) { weight = pTSWeights->GetSampledValue(nn); }
      if((obsval==RAV_BLANK_DATA) || (modval==RAV_BLANK_DATA)) { weight=0.0; }

      sum1 += weight*fabs(obsval - modval);
      sum2 += weight*fabs(avgobs - modval);
      N   += weight;
    }
    if (N>0.0)
    {
      return (sum1 / sum2);
    }
    else
    {
      string warn = "DIAG_RABSERR not calculated. Missing non-zero weighted observations during simulation duration.";
      WriteWarning(warn, Options.noisy);
      return -ALMOST_INF;
    }
    break;
  }
  case(DIAG_PERSINDEX)://----------------------------------------------------
  {
    bool ValidObs = false;
    bool ValidMod = false;
    bool ValidWeight = false;

    double sum1(0.0),sum2(0.0);
    N=0;
    double N1(0.0),N2(0.0);

    for(nn=nnstart+1; nn<nnend; nn++)
    {
      double prvmodval = pTSMod->GetSampledValue(nn - 1);
      double prvweight=1;
      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      weight=1.0;
      if (pTSWeights != NULL) { weight = pTSWeights->GetSampledValue(nn); prvweight = pTSWeights->GetSampledValue(nn - 1); }

      if (obsval != RAV_BLANK_DATA) { ValidObs = true; }
      if (modval != RAV_BLANK_DATA && prvmodval != RAV_BLANK_DATA) { ValidMod = true; }
      if (weight != 0 && prvweight != 0) { ValidWeight = true; }

      if (obsval != RAV_BLANK_DATA && modval != RAV_BLANK_DATA && prvmodval != RAV_BLANK_DATA && weight != 0 && prvweight != 0)
      {
        sum1 +=weight*(pow( modval-obsval, 2));
        sum2 +=prvweight*(pow( modval-prvmodval, 2));
        ++N;
        // consider alternative weighting N scheme
        // sum1 +=weight*(pow( modval-obsval, 2);
        // sum2 +=prvweight*(pow( modval-prvmodval, 2);
        // N1 += weight
        // N2 += 0.5*(prvweight+weight)
      }
    }
    
    if(N>0.0) {
      if (ValidObs && ValidMod && ValidWeight)
      {
        return 1 - (sum1 / sum2);
        // return 1 - ( (sum1/N1) / (sum2/N2) )
      }
      else
      {
        string p1 = "OBSERVED_DATA  "+filename;
        string p2 = "MODELED_DATA  ";
        string p3 = "WEIGHTS ";
        if (ValidObs) { p1 = ""; }
        if (ValidMod) { p2 = ""; }
        if (ValidWeight) { p3 = ""; }

        string warn = "DIAG_PERSINDEX not performed correctly. Check " + p1 + p2 + p3;
        WriteWarning(warn, Options.noisy);
        return -ALMOST_INF;
      }
    break;
    }
    else
    {
      string warn = "DIAG_PERSINDEX not not calculated. Missing non-zero weighted observations during simulation duration.";
      WriteWarning(warn,Options.noisy);
      return -ALMOST_INF;
    }
    break;
  }
  case(DIAG_NSE4)://----------------------------------------------------
  {
    double avgobs=0.0;
    N=0;

    for(nn=nnstart;nn<nnend;nn++)
    {
      weight=1.0;
      obsval = pTSObs->GetSampledValue(nn);
      if(pTSWeights != NULL) { weight=pTSWeights->GetSampledValue(nn); }
      if(obsval==RAV_BLANK_DATA) { weight=0.0; }
      avgobs+=weight*obsval;
      N     +=weight;
    }
    if(N>0.0) { avgobs/=N; }

    double sum1(0.0),sum2(0.0);
    for(nn=nnstart;nn<nnend;nn++)
    {
      weight=1.0;
      obsval = pTSObs->GetSampledValue(nn);
      modval = pTSMod->GetSampledValue(nn);
      if(pTSWeights != NULL) { weight=pTSWeights->GetSampledValue(nn); }
      if((obsval==RAV_BLANK_DATA) || (modval==RAV_BLANK_DATA)) { weight=0.0; }

      sum1 += weight*pow(obsval - modval,4);
      sum2 += weight*pow(obsval - avgobs,4);
    }

    if(N>0)
    {
      return 1.0 - (sum1 / sum2);
    }
    else
    {
      string warn = "DIAG_NSE4 not calculated. Missing non-zero weighted observations during simulation duration.";
      WriteWarning(warn,Options.noisy);
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
/*****************************************************************
Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/
//////////////////////////////////////////////////////////////////
/// \brief Implementation of the CDiagnosticPeriod constructor/destructor
/// \param name [in] period name (e.g., "CALIBRATION")
/// \param startdate [in] string date of period beginning "yyyy-mm-dd" format
/// \param enddate [in] string date of period ending "yyyy-mm-dd" format
/// \param Options [in] options structure
//
CDiagPeriod::CDiagPeriod(string name,string startdate,string enddate,const optStruct &Options)
{
  _name=name;
  time_struct tt;
  tt=DateStringToTimeStruct(startdate,"00:00:00",Options.calendar);
  _t_start = TimeDifference(Options.julian_start_day,Options.julian_start_year,tt.julian_day,tt.year,Options.calendar);
  
  tt=DateStringToTimeStruct(enddate  ,"00:00:00",Options.calendar);
  _t_end=   TimeDifference(Options.julian_start_day,Options.julian_start_year,tt.julian_day,tt.year,Options.calendar);

  if (_t_end<=_t_start){WriteWarning("CDiagPeriod: :EvaluationPeriod startdate after enddate.",Options.noisy); }
}
CDiagPeriod::~CDiagPeriod() {}
//////////////////////////////////////////////////////////////////
/// \brief Implementation of the CDiagnosticPeriod accessors
//
string   CDiagPeriod::GetName()      const{return _name;}
double   CDiagPeriod::GetStartTime() const{return _t_start;}
double   CDiagPeriod::GetEndTime()   const{return _t_end;}
