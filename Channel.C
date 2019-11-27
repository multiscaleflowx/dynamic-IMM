#include "Channel.H"
#include "common.H"

namespace cfdsim {

  Channel::Channel() {}

  Channel::Channel(double len, double heightAtEnd, double heightAtNeck) {
    l = len;
    hEnd = heightAtEnd;
    hNeck = heightAtNeck;
  }

  // The channel shrinks linearly until it reaches the neck
  // and then expands to its original size.
  // s is normalised.
  double Channel::heightAt(double s) {
    if(s < 0.0) {
      haltMPMD(" the value supplied to heightAt was negative.");
    }
    if(s >= 1.0) {
      haltMPMD("the value supplied to heightAt was one or greater.");
    }
    if(s <= 0.5) {
      return (1 - 2*s)*hEnd + 2*s*hNeck;
    }
    else {
      return (2 - 2*s)*hNeck + (2*s - 1)*hEnd;
    }
  }

  double Channel::length() {
    return l;
  }

}
