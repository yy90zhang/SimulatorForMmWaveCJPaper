#ifndef NETWORK_H_INCLUDED
#define NETWORK_H_INCLUDED

#include "randGenerator.h"
#define NOISE pow(10,-11.4)//0.0001
#define ALPHA_N 4        ///path loss exponent of NLoS links
#define ALPHA_L 2        ///path loss exponent of LoS links
#define P_LOS 0.2      ///LoS probability
#define R_LOS 200        ///radius of the LoS ball
#define R0 100           ///D2D link length <= R_LOS
#define NL 3
#define NN 2
#define PI 3.1415926
#define TRANSMITTER 0
#define RECEIVER 1
#define JAMMER 2
#define EAVESDROPPER 3
using namespace std;

/*-------------------struct of node's position------------------*/
typedef struct _position{
    double r;                       ///coordinate x
    double theta;                       ///coordinate y
}PT;

typedef struct _antenna{
    double AOB;                    ///angle of boresight
    double beamwidth_m;            ///beamwidth of the main lobe
    double back_lobe_gain;         ///back lobe gain
    double main_lobe_gain;         ///main lobe gain
}AT;

typedef struct _pair{   ///association between a Tx and Rx, used for identifying the associated active rx of a non-active tx
    int tx;             ///index of the transmitter
    int rx;             ///index of the receiver
}PAIR;

/**get the random channel gain and path loss of a link of distance d**/
void randomPathLossAndChannelGain(double d, double & PL, double & gain);

/**vector addition in polar coordinate system, used for:
(1)determining the position of the d2d receiver of a transmitter
(2)shifting a node with respect to a reference node
(i.e., changing the origin of the coordinate system)**/
/*pt: point to shift; ref_pt: reference point; returned value: pt1 + pt2*/
PT vectorAddition(PT pt1, PT pt2);


/*--------------------Class Node-----------------------------*/
class Node{
private:
    PT position;                        ///node position
    double power;                       ///transmit power
    double decodeTh;                    ///SINR threshold for successful decoding
    double SINR;                        ///received SINR
    int dFlag;                          ///flag indicating a successful decoding by 1, and a failed decoding by 0
    int index;                          ///index of node
    AT antenna;                             ///antenna


public:
    Node();                             ///constructor
    void setPosition(PT ps);            ///set node's position
    void setPower(double tp);           ///set transmit power
    void setDecodeTH(double dt);        ///set decoding threshold
    void setSINR(double rcvP, double inter, double noise); ///calculate the received SINR based on received power (rcvP), interference (inter) and noise (noise)
    void setAntenna(AT a);
    double transmit();                  ///transmit message, in fact only power is of interest
    void decode();                      ///decode message
    PT getPosition();                   ///get node's position
    AT getAntenna();
    double dis2RefNode(Node ref_node);  ///return its distance to a reference node
    int getDecodeFlag();
    double getSINR();
    double getDecodeTH();
    ~Node();                            ///deconstructor
};

/*--------------class Eavesdropper inherited from class Node-------*/
class Eavesdropper:public Node{
};


/*-------------class Network-------------------------------*/
class Network{
private:
    /**network parameters**/
    double radius;                             ///radius of network disk int num_D;
    int num_D;                                 ///number of D2D pairs (Poisson random variable with mean mean_D)
    int num_PJ;                                ///number of potential jammers
    int num_J;                                 ///number of jammers
    int num_E;                                 ///number of eavesdroppers
    double mean_D;                             ///expected number of legitimate nodes in the network
    double mean_E;                             ///expected number of eavesdroppers in the network
    double mean_PJ;                            ///expected number of potential jammers
    double mean_J;                             ///expected number of jammers
    double density_D;                          ///density of legitimate nodes (lambda_T)
    double density_E;                          ///density of eavesdroppers (lambda_E)
    double density_PJ;                         ///density of potential jammers (lambda_J)
    Node * transmitters;                       ///set of D2D transmitters
    Node * receivers;                          ///set of D2D receivers (tx-rx pairing: 0<->0, 1<->1, 2<->2......)
    Node * potentialJammers;                   ///set of potentail jammers
    int * index_Jammers;                       ///index of jammers
    int * associatedReceivers;                 ///index of associated receiver of each jammer
    Eavesdropper *eavesdroppers;               ///set of eavesdroppers
    /**LDCJ parameters**/
    double q;                                  ///probability of jammer selection probability in LoS balls

public:
   Network();     ///constructor with default parameters
   void setRadius(double r);                                                                             ///set network size
   void setDensity_D(double den_D);                                                                                ///set density of legitimate nodes and eavesdroppers
   void setDensity_E(double den_E);
   void setDensity_PJ(double den_PJ);
   void set_q(double q);
   void initD2DPairs(double tp, double w_t, double g_b_t, double g_m_t, double w_r, double g_b_r, double g_m_r, double dt, int flag);
   void initJammers(double tp, double w_t, double g_b_t, double g_m_t);
   void initEvesdroppers(double w_e, double g_b_e, double g_m_e, double dt);
   void deleteD2DPairs();
   void deleteJammers();
   void deleteEavedroppers();
   double interference2RefNode(int index, int type_interfer, int type_receiver);
   bool legitmateTransmission();
   bool eavesdropping();
   int getNumD();
};


/*-----------------------struct paramters of network--------------------*/
typedef struct _params{
    /*network-related parameters*/
    double NW_R;                        ///network radius
    double den_d;                       ///density of D2D pair
    double den_pj;                      ///density of potential jammers
    double den_e;                       ///density of eavesdroppers
    /*transmission-related parameters*/
    double tp;                          ///D2D transmit power
    double dt_e;                        ///eavesdropper decoding threshold
    double dt_d;                        ///legitimate node decoding threshold
    /*LDCJ parameters*/
    double q;                           ///probability of jammer selection probability in LoS balls
    /*antenna-related paramters*/
    double w_t;                         ///tx's beamwidth of main lobe
    double w_r;                         ///rx's beamwidth of main lobe
    double w_e;                         ///eve's beamwidth of main lobe
    double g_b_t;                       ///tx's back lobe gain
    double g_b_r;                       ///rx's back lobe gain
    double g_b_e;                       ///eve's back lobe gain
    double g_m_t;                       ///tx's main lobe gain
    double g_m_r;                       ///rx's main lobe gain
    double g_m_e;                       ///eve's main lobe gain
}Parameters;

/*-----------------------class Simulator-------------------------------*/
class Simulator{
private:
    Network netW;
    int loop_num;
    Parameters pmts;

public:
    Simulator();
    void setLoopNum(int ln);
    void initParameters(Parameters pm);
    void creatNWInstance();
    void LTInterferenceOfTypicalD2DReceiver(double s);            ///Laplace Transform of total interference at the typical D2D receiver
    double connectionProbability();
    double secrecyProbabiltiy();
    ~Simulator();

};
/**check if node from is in the main lobe of node to**/
bool isInMainLobe(Node from, Node to);
double randomAntennaGain(Node tx, Node rx);
#endif // NETWORK_H_INCLUDED
