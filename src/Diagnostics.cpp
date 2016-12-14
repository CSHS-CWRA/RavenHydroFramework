/*----------------------------------------------------------------
Raven Library Source Code
Copyright © 2008-2015 the Raven Development Team, Ayman Khedr, Konhee Lee
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
    case(DIAG_NASH_SUTCLIFFE):{return "DIAG_NASH_SUTCLIFFE"; break;}
    case(DIAG_RMSE):          {return "DIAG_RMSE";break;}
    case(DIAG_PCT_BIAS):      {return "DIAG_PCT_BIAS";break;}
    case(DIAG_ABSERR):        {return "DIAG_ABSERR";break;}
    case(DIAG_ABSMAX):        {return "DIAG_ABSMAX";break;}
    case(DIAG_PDIFF):         {return "DIAG_PDIFF";break;}
    case(DIAG_TMVOL):         {return "DIAG_TMVOL";break;}
    case(DIAG_RCOEF):         {return "DIAG_RCOEF"; break;}
    case(DIAG_NSC):           {return "DIAG_NSC";break;}
    case(DIAG_RSR):           {return "DIAG_RSR";break;}
    case(DIAG_R2):            {return "DIAG_R2";break;}
    case(DIAG_CUMUL_FLOW):    {return "DIAG_CUMUL_FLOW";break;}
    case(DIAG_LOG_NASH):      {return "DIAG_LOG_NASH";break;}
    case(DIAG_KLING_GUPTA):   {return "DIAG_KLING_GUPTA";break;}
    default:                  {return "";break;}
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
  int nVals=pTSObs->GetNumSampledValues();
  double obsval,modval;
  double weight=1;
  int skip=0;
  if (!strcmp(pTSObs->GetName().c_str(), "HYDROGRAPH") && (Options.ave_hydrograph == true)){ skip = 1; }//skips first point (initial conditions)

  switch (_type)
  {
    case(DIAG_NASH_SUTCLIFFE)://-----------------------------------------
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

        if (obsval != CTimeSeriesABC::BLANK_DATA) { ValidObs = true; }
        if (modval != CTimeSeriesABC::BLANK_DATA) { ValidMod = true; }
        if (weight != 0) { ValidWeight = true; }
        
        if (obsval != CTimeSeriesABC::BLANK_DATA && modval != CTimeSeriesABC::BLANK_DATA && weight != 0) {
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

          if (obsval != CTimeSeriesABC::BLANK_DATA) { ValidObs = true; }
          if (modval != CTimeSeriesABC::BLANK_DATA) { ValidMod = true; }
          if (weight != 0) { ValidWeight = true; }

          if (obsval != CTimeSeriesABC::BLANK_DATA && modval != CTimeSeriesABC::BLANK_DATA && weight != 0) {
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
        string p1 = "OBSERVED_DATA  ";
        string p2 = "MODELED_DATA  ";
        string p3 = "WEIGHTS ";
        if (ValidObs) { p1 = ""; }
        if (ValidMod) { p2 = ""; }
        if (ValidWeight) { p3 = ""; }

        string warn = "DIAG_NASH_SUTCLIFFE not performed correctly. Check " + p1 + p2 + p3;
        WriteWarning(warn, Options.noisy);
        return -ALMOST_INF;
      }
      break;
    }
    case(DIAG_RMSE)://----------------------------------------------------
    {
        bool ValidObs = false;
        bool ValidMod = false;
        bool ValidWeight = false;

        double sum;
        N=0;
        sum=0;

        for (nn=skip;nn<nVals;nn++)
        {
            obsval=pTSObs->GetSampledValue(nn);
            modval=pTSMod->GetSampledValue(nn);
            weight=1.0;
            if(pTSWeights != NULL){ weight=pTSWeights->GetSampledValue(nn);}

            if (obsval != CTimeSeriesABC::BLANK_DATA) { ValidObs = true; }
            if (modval != CTimeSeriesABC::BLANK_DATA) { ValidMod = true; }
            if (weight != 0) { ValidWeight = true; }

            if (obsval != CTimeSeriesABC::BLANK_DATA && modval != CTimeSeriesABC::BLANK_DATA && weight != 0) {
                sum+=pow(obsval-modval,2)*weight;
                N+=weight;
            }
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

            string warn = "DIAG_RMSE not performed correctly. Check " + p1 + p2 + p3;
            WriteWarning(warn, Options.noisy);
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

            if (obsval != CTimeSeriesABC::BLANK_DATA) { ValidObs = true; }
            if (modval != CTimeSeriesABC::BLANK_DATA) { ValidMod = true; }
            if (weight != 0) { ValidWeight = true; }

            if (obsval != CTimeSeriesABC::BLANK_DATA && modval != CTimeSeriesABC::BLANK_DATA && weight != 0) {
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
            string p1 = "OBSERVED_DATA  ";
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

            if (obsval != CTimeSeriesABC::BLANK_DATA) { ValidObs = true; }
            if (modval != CTimeSeriesABC::BLANK_DATA) { ValidMod = true; }
            if (weight != 0) { ValidWeight = true; }

            if (obsval != CTimeSeriesABC::BLANK_DATA && modval != CTimeSeriesABC::BLANK_DATA && weight != 0){
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
            string p1 = "OBSERVED_DATA  ";
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
            
            if (obsval != CTimeSeriesABC::BLANK_DATA) { ValidObs = true; }
            if (modval != CTimeSeriesABC::BLANK_DATA) { ValidMod = true; }
            if (weight != 0) { ValidWeight = true; }

            if (obsval != CTimeSeriesABC::BLANK_DATA && modval != CTimeSeriesABC::BLANK_DATA && weight != 0) {
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
            string p1 = "OBSERVED_DATA  ";
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

            if (obsval != CTimeSeriesABC::BLANK_DATA) { ValidObs = true; }
            if (modval != CTimeSeriesABC::BLANK_DATA) { ValidMod = true; }
            if (weight != 0) { ValidWeight = true; }

            if (obsval != CTimeSeriesABC::BLANK_DATA && modval != CTimeSeriesABC::BLANK_DATA && weight != 0) {
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
            string p1 = "OBSERVED_DATA  ";
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

            if (obsval != CTimeSeriesABC::BLANK_DATA && modval != CTimeSeriesABC::BLANK_DATA && weight != 0)
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

            if (obsval != CTimeSeriesABC::BLANK_DATA) { ValidObs = true; }
            if (modval != CTimeSeriesABC::BLANK_DATA) { ValidMod = true; }
            if (weight != 0) { ValidWeight = true; }

            if (obsval != CTimeSeriesABC::BLANK_DATA && modval != CTimeSeriesABC::BLANK_DATA && weight != 0)
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
            string p1 = "OBSERVED_DATA  ";
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

            if (obsval != CTimeSeriesABC::BLANK_DATA && nxtobsval != CTimeSeriesABC::BLANK_DATA) { ValidObs = true; }
            if (modval != CTimeSeriesABC::BLANK_DATA && nxtmodval != CTimeSeriesABC::BLANK_DATA) { ValidMod = true; }
            if (weight != 0 && nxtweight != 0) { ValidWeight = true; }

            if (obsval != CTimeSeriesABC::BLANK_DATA && modval != CTimeSeriesABC::BLANK_DATA && nxtobsval != CTimeSeriesABC::BLANK_DATA && nxtmodval != CTimeSeriesABC::BLANK_DATA && weight != 0 && nxtweight != 0)
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

            if (obsval != CTimeSeriesABC::BLANK_DATA && nxtobsval != CTimeSeriesABC::BLANK_DATA) { ValidObs = true; }
            if (modval != CTimeSeriesABC::BLANK_DATA && nxtmodval != CTimeSeriesABC::BLANK_DATA) { ValidMod = true; }
            if (weight != 0 && nxtweight != 0) { ValidWeight = true; }

            if (obsval != CTimeSeriesABC::BLANK_DATA && modval != CTimeSeriesABC::BLANK_DATA && nxtobsval != CTimeSeriesABC::BLANK_DATA && nxtmodval != CTimeSeriesABC::BLANK_DATA && weight != 0 && nxtweight != 0)
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
            string p1 = "OBSERVED_DATA  ";
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

            if (obsval != CTimeSeriesABC::BLANK_DATA && nxtobsval != CTimeSeriesABC::BLANK_DATA) { ValidObs = true; }
            if (modval != CTimeSeriesABC::BLANK_DATA && nxtmodval != CTimeSeriesABC::BLANK_DATA) { ValidMod = true; }
            if (weight != 0 && nxtweight != 0) { ValidWeight = true; }

            if (obsval != CTimeSeriesABC::BLANK_DATA && modval != CTimeSeriesABC::BLANK_DATA && nxtobsval != CTimeSeriesABC::BLANK_DATA && nxtmodval != CTimeSeriesABC::BLANK_DATA && weight != 0 && nxtweight != 0)
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
            string p1 = "OBSERVED_DATA  ";
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

            if (obsval != CTimeSeriesABC::BLANK_DATA) { ValidObs = true; }
            if (modval != CTimeSeriesABC::BLANK_DATA) { ValidMod = true; }
            if (weight != 0) { ValidWeight = true; }

            if (obsval != CTimeSeriesABC::BLANK_DATA && modval != CTimeSeriesABC::BLANK_DATA && weight != 0)
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

            if (obsval != CTimeSeriesABC::BLANK_DATA) { ValidObs = true; }
            if (modval != CTimeSeriesABC::BLANK_DATA) { ValidMod = true; }
            if (weight != 0) { ValidWeight = true; }

            if (obsval != CTimeSeriesABC::BLANK_DATA && modval != CTimeSeriesABC::BLANK_DATA && weight != 0)
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
            string p1 = "OBSERVED_DATA  ";
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

            if (obsval != CTimeSeriesABC::BLANK_DATA) { ValidObs = true; }
            if (modval != CTimeSeriesABC::BLANK_DATA) { ValidMod = true; }
            if (weight != 0) { ValidWeight = true; }

            if (obsval != CTimeSeriesABC::BLANK_DATA && modval != CTimeSeriesABC::BLANK_DATA && weight != 0)
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

            if (obsval != CTimeSeriesABC::BLANK_DATA) { ValidObs = true; }
            if (modval != CTimeSeriesABC::BLANK_DATA) { ValidMod = true; }
            if (weight != 0) { ValidWeight = true; }

            if (obsval != CTimeSeriesABC::BLANK_DATA && modval != CTimeSeriesABC::BLANK_DATA && weight != 0)
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
            string p1 = "OBSERVED_DATA  ";
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
        if (obsval <= 0.0){obsval=CTimeSeriesABC::BLANK_DATA;} //negative values treated as invalid observations/modeled
        else              {obsval=log(obsval);               }
        if (modval <= 0.0){modval=CTimeSeriesABC::BLANK_DATA;}
        else              {modval=log(modval);               }

        if (obsval != CTimeSeriesABC::BLANK_DATA) { ValidObs = true; }
        if (modval != CTimeSeriesABC::BLANK_DATA) { ValidMod = true; }
        if (weight != 0) { ValidWeight = true; }
        
        if ((obsval != CTimeSeriesABC::BLANK_DATA) && (modval != CTimeSeriesABC::BLANK_DATA) && (weight != 0)) 
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
        if (obsval <= 0.0){obsval=CTimeSeriesABC::BLANK_DATA;} //negative values treated as invalid observations/modeled
        else              {obsval=log(obsval);               }
        if (modval <= 0.0){modval=CTimeSeriesABC::BLANK_DATA;}
        else              {modval=log(modval);               }

        if ((obsval != CTimeSeriesABC::BLANK_DATA) && (modval != CTimeSeriesABC::BLANK_DATA) && (weight != 0)) {
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
        string p1 = "OBSERVED_DATA  ";
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

        if (obsval != CTimeSeriesABC::BLANK_DATA) { ValidObs = true; }
        if (modval != CTimeSeriesABC::BLANK_DATA) { ValidMod = true; }
        if (weight != 0) { ValidWeight = true; }
        
        if ((obsval != CTimeSeriesABC::BLANK_DATA) && (modval != CTimeSeriesABC::BLANK_DATA) && (weight != 0)) 
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
        string p1 = "OBSERVED_DATA  ";if (ValidObs) { p1 = ""; }
        string p2 = "MODELED_DATA  "; if (ValidMod) { p2 = ""; }
        string p3 = "WEIGHTS ";       if (ValidWeight) { p3 = ""; }

        string warn = "DIAG_CUMUL_FLOW not performed correctly. Check " + p1 + p2 + p3;
        WriteWarning(warn, Options.noisy);
        return ALMOST_INF;
      }
      break;
    }
  default:
    {
      return 0.0; break;
    }
  }//end switch


}