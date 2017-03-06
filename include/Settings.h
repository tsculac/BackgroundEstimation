#ifndef Settings_h
#define Settings_h

using namespace std;

class Settings
{

public:
	
	Settings();
	~Settings();
   
   enum _process { Data = 0, WZ = 1, MAX_NUM_OF_PROCESSES };
   enum _flavour { ele = 0, mu = 1, MAX_NUM_OF_FLAVOURS };
   enum _final_state { fs4mu = 0, fs4e = 1, fs2e2mu = 2, fs2mu2e = 3, fs4l = 4, MAX_NUM_OF_FINAL_STATES };
   enum _eta_bins { EB = 0, EE = 1, MAX_NUM_OF_ETA_BINS};
   
   static const int num_of_processes         = MAX_NUM_OF_PROCESSES;
   static const int num_of_flavours          = MAX_NUM_OF_FLAVOURS;
   static const int num_of_final_states      = MAX_NUM_OF_FINAL_STATES;
   static const int num_of_eta_bins          = MAX_NUM_OF_ETA_BINS;

   private:
      
};
#endif
