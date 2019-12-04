#include <map>
#include <string>
#include <istream>
#include <ostream>
#include <fstream>
#include <set>          // std::set
#include <vector>       // std::vector

#include "CFDSim.H"
#include "common.H"

namespace cfdsim {

  void CFDSim::initialise(int nMicro_) {
    Info << nl << "Initialising macro solver" << nl << endl;

    double F_ = initialValues["force"];

    s_.setSize(nMicro_, 0.0);
    mDot_.setSize(nMicro_, 0.0);
    k_.setSize(nMicro_, 1.0);
    f_.setSize(nMicro_, F_); // F_ is the correct initial value for f_.
    f_old_.setSize(nMicro_, 0.0);
    phi_.setSize(nMicro_, 0.0);
    phiCoeffs_.setSize(nMicro_, 0.0);
    kIs_.setSize(nMicro_);

    std::set<double> normalisedPositions;
    std::cout << "normalisedPositions:";
    for(Region r : regions) {
      normalisedPositions.insert(r.sNorm);
      std::cout << " " << r.sNorm;
    }
    std::cout << std::endl;

    // Compute the heights of the channel for the various regions.
    label i = 0;
    for(double s : normalisedPositions) {
      s_[i] = s*channel.length();
      std::cout << "s_[" << i << "] = " << s_[i] << std::endl;
      i++;
    }

    forAll(kIs_, i) {
      kIs_[i].setSize(nIter_);
      forAll(kIs_[i], j) {
	kIs_[i][j].setSize(3, 0.0);
      }
    }
    Info << nl << "Macro solver initialised" << nl << endl;
  }

  void CFDSim::solve(label iter_, label nMicro_) {
    scalar pi = constant::mathematical::pi;
    double F_ = initialValues["force"];

    Info << nl << "iter:" << iter_
	 << " solving macro equation"
	 << endl;

    // 1. set k_i's from current MD simulation data

    if(iter_ == 0) {
      Info << nl << "applying initial flow resistance adjustments" << endl;

      for(label i = 0; i < nMicro_; i++) {
	ofs << "mDot_[" << i << "] = " << mDot_[i] << std::endl;
	if(mDot_[i] != 0.0) {
	  k_[i] = f_[i]/mDot_[i];

	  kIs_[i][iter_][0] = mDot_[i];
	  kIs_[i][iter_][1] = f_[i];
	  kIs_[i][iter_][2] = 0.1*mDot_[i]; //mDotError_[i] //*** need to pass flow rate error
	}
	ofs << "kIs_[" << i << "][" << iter_ << "][0] = " << kIs_[i][iter_][0] << std::endl;
	ofs << "kIs_[" << i << "][" << iter_ << "][1] = " << kIs_[i][iter_][1] << std::endl;
	ofs << "kIs_[" << i << "][" << iter_ << "][2] = " << kIs_[i][iter_][2] << std::endl;
	ofs << "k_[" << i << "] = " << k_[i] << std::endl;
      }
    }
    else {
      Info << nl << "Calculating k_i's" << nl << endl;
 
      for(label i = 0; i < nMicro_; i++) {
	kIs_[i][iter_][0] = mDot_[i];
	kIs_[i][iter_][1] = f_[i];
	kIs_[i][iter_][2] = 0.1*mDot_[i]; //mDotError_[i] //*** need to pass flow rate error

	// find k_i's
    
	scalar sumTop = 0.0;
	scalar sumBottom = 0.0;

	for (label j=0; j < nIter_; j++) {
	  if(j < iter_) {
	    scalar mDot = kIs_[i][j][0];
	    scalar force = kIs_[i][j][1];
	    scalar mDotError = kIs_[i][j][2];

	    sumTop += force*force/(mDotError*mDotError);
	    sumBottom += force*mDot/(mDotError*mDotError);
	  }
	}

	k_[i] = sumTop/sumBottom;

	ofs << "kIs_[" << i << "][" << iter_ << "][0] = " << kIs_[i][iter_][0] << std::endl;
	ofs << "kIs_[" << i << "][" << iter_ << "][1] = " << kIs_[i][iter_][1] << std::endl;
	ofs << "kIs_[" << i << "][" << iter_ << "][2] = " << kIs_[i][iter_][2] << std::endl;
	ofs << "k_[" << i << "] = " << k_[i] << std::endl;
	//Calculate ki errors here //****
      }
    }

    // 2. set LU matrix (this is the hardest and trickiest part)

    std::cout << "Creating matrix" << std::endl;
    simpleMatrix<scalar> luMatrix(nMicro_+1, 0.0, 0.0);

    // column configuration

    // phiCoeff's... mbar

    label r = 0; // row counter

    label c = 0; // column pusher

    // overall pressure drop integral equation

    double L_ = channel.length();

    luMatrix[r][c+0] = L_;

    for(label j=1; j <= (nMicro_-1)/2; j++) {
      luMatrix[r][c+(2*j)-1] = (L_/(2*pi*j))*(1-Foam::cos(2*pi*j));
      luMatrix[r][c+(2*j)+1-1] = L_*(Foam::sin(2*pi*j))/(2*pi*j);
    }

    luMatrix.source()[r]=0.0;

    r++;

    // flow response equations

    for(label i = 0; i < nMicro_; i++) {
      c = nMicro_;

      luMatrix[r][c]=-rhoN_*k_[i];

      c = 0;

      luMatrix[r][c+0]=1;

      for(label j=1; j <= (nMicro_-1)/2; j++) {
	luMatrix[r][c+(2*j)-1] = Foam::sin(2*pi*s_[i]*j/L_);
	luMatrix[r][c+(2*j)+1-1] = Foam::cos(2*pi*s_[i]*j/L_);
      }

      luMatrix.source()[r]=-rhoN_*F_;

      r++;
    }

    //... done

    // 3. solve matrix
    scalarField M = luMatrix.LUsolve();

    // 4. set new phi's

    for(label i = 0; i < nMicro_; i++) {
      phiCoeffs_[i] = M[i];
    }

    for(label i = 0; i < nMicro_; i++) {
      phi_[i] = phiCoeffs_[0];

      for(label j=1; j <= (nMicro_-1)/2; j++) {
	phi_[i] += phiCoeffs_[(2*j)-1]*Foam::sin(2*pi*s_[i]*j/L_);
	phi_[i] += phiCoeffs_[(2*j)+1-1]*Foam::cos(2*pi*s_[i]*j/L_);
      }
    }

    //5. Set new forces

    List<scalar> mDotMacro_(nMicro_,0.0); //predicted Mass Flow Rates

    std::cout << "Setting new forces" << std::endl;
    for(label i = 0; i < nMicro_; i++) {
      f_old_[i] = f_[i];
      f_[i] = (F_ + phi_[i]/rhoN_);
      ofs << "f_[" << i << "] = " << f_[i] << " = " << F_ << " + " << phi_[i] << "/" << rhoN_ << std::endl;

      if(k_[i] != 0.0) {
	mDotMacro_[i] = f_[i]/k_[i];
      }
    }

    // next iteration

    //iter_++; // important

    Info << nl << "Statistics" << endl;
    Info << nl << "id \t s_i \t f_i \t mDotMicro_i \t mDotMacro_i \t k_i"
	 << "\t phi_i \t phiCoeffs_i \t"
	 << endl;

    for(label i = 0; i < nMicro_; i++) {
      Info << i << "\t" << s_[i] << "\t" << f_[i] << "\t" << mDot_[i]
	   << "\t" << mDotMacro_[i] << "\t" << k_[i] << "\t" << phi_[i]
	   << "\t" << phiCoeffs_[i]
	   << endl;
    }

    ofs << std::endl << "Statistics" << std::endl;
    ofs << std::endl << "id \t s_i \t f_i \t mDotMicro_i \t mDotMacro_i \t k_i"
        << "\t phi_i \t phiCoeffs_i \t"
        << std::endl;

    for(label i = 0; i < nMicro_; i++) {
      ofs << i << "\t" << s_[i] << "\t" << f_[i] << "\t" << mDot_[i]
	  << "\t" << mDotMacro_[i] << "\t" << k_[i] << "\t" << phi_[i]
	  << "\t" << phiCoeffs_[i]
	  << std::endl;
    }
  }

  void CFDSim::calculateOutputs(int64_t timestep) {
    int32_t posn = 0;
    for(std::string ifn : interfaceNames) {
      for(std::string v : outVars) {
	if(v == std::string("force")) {
	  std::cout << "f_[" << posn << "] = " << f_[posn] << std::endl;
	  std::cout << "f_old_[" << posn << "] = " << f_old_[posn] << std::endl;
	  double deltaF_ = f_[posn];
	  outVarValues[ifn][v] = deltaF_;
	  std::cout << "Force delta set to " << deltaF_ << std::endl;
	}
	else if(v == std::string("cfd_push_param")) {
	  outVarValues[ifn][v] = double(timestep);
	}
	else {
	  haltMPMD("unexpected variable.");
	}
      }
      posn++;
    }
  }

  void CFDSim::setInputs() {
    // This relies on the regions being correctly ordered.
    sort(regions.begin(), regions.end());
    int32_t j = 0;
    for(Region r : regions) {
      std::cout << "Sorted: ifn = " << r.interfaceName << std::endl;
      mDot_[j] = averages[r.interfaceName][std::string("mass_flow_x")].back();
      std::cout << "equib: mDot_[" << j << "] = " << mDot_[j] << std::endl;
      j++;
    }
  }

}
