#include <assert.h>

#include "Channel.H"

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
    assert(s >= 0.0);
    assert(s < 1.0);
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
