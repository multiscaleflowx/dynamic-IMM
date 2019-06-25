#ifndef CHANNEL_H
#define CHANNEL_H

namespace cfdsim {
  class Channel {
    double l;
    double hEnd;
    double hNeck;

  public:
    Channel();

    Channel(double, double, double);

    double heightAt(double);

    double length();
  };
}

#endif
