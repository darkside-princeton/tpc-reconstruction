/*---------------------------------------------------------------------*
| Event.h                                                              |
|                                                                      |
| Description: A class that holds an vector of channel objects and     |
|              performs analysis on clusters detected across channels. |   
|              Sets up most of the output variables/global variables   |
|              held in the Output class.                               |
|                                                                      |
*----------------------------------------------------------------------*/

#ifndef Event_h
#define Event_h

#include <iostream>
#include <map>
#include "Config.h"
#include "Pulse.h"
#include "Channel.h"
#include "Cluster.h"

using namespace std;

class Event{
 public:
  Event() {
    Event(-1);
  }

  Event(int ev_id){
    _eventID = ev_id;
    for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
      Channel* ch = new Channel(Config::Get()->GetParameterIvec("ch_ids")[i_ch]);
      _channels.push_back(ch);
    }
    _sum_wf = NULL;
    _cluster_overlap = false;
  }

  ~Event() { ; }

  // Get Functions
  int              GetEventID()        {return _eventID; }
  Channel*         GetChannel(int c)   {return _channels[c];}
  vector<Channel*> GetChannels()       {return _channels;}
  vector<Cluster*> GetClusters()       {return _clusters;}
  Waveform*        GetSumWF()          {return _sum_wf;}
  bool             GetClusterOverlap() {return _cluster_overlap;}


  // Set Functions
  void SetRawWF(int channel, Waveform*  wf);
  void SetEventID(int ev_id)                   {_eventID = ev_id;}
  void SetClusterOverlap(bool cluster_overlap) {_cluster_overlap = cluster_overlap;}

  // Processing Functions
  void ProcessEvent();
  void ProcessEventNoise();
  void FindClusters();
  void FindClusterPulses();
  void ClearClusters() {_clusters.clear();}
  void ClusterDeconv(Waveform* response);
  void CalcSumWaveform();


  // Output Functions
  void PulseOutput();
  void ClusterOutput();

  // Display Functions
  void EventPrint(int channel);
  void DisplayEvent(int channel, double min_x = 0, double max_x = Config::Get()->GetParameterI("num_samps") * 1e-3 / Config::Get()->GetParameterD("sampling_rate"), 
		    double min_y = 0, double max_y = 0);
  
 
  void EventFree();

private:
  int              _eventID;
  Waveform*        _sum_wf;  
  vector<Cluster*> _clusters;
  bool             _cluster_overlap;
  vector<Channel*> _channels;


  // Display Functions
  void PlotClusterTimes(int channel, double min_x, double max_x, double min_y, double max_y);
  void PlotClusterPulses(int channel, double min_x, double max_x, double min_y, double max_y);
  void PlotSumWaveform(double min_x, double max_x, double min_y, double max_y);

};



#endif
