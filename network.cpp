#include "network.h"
#include <iostream>
#include <math.h>

void randomPathLossAndChannelGain(double d, double & pl, double & gain)
{
    double alpha = ALPHA_N;
    double shape = NN, scale = (double)1/NN;
    if (d <= R_LOS && uniform_rv_generator(0.0, 1.0) < P_LOS)
    {//if the link is LOS
        alpha = ALPHA_L;
        shape = NL;
        scale = (double)1/NL;
    }
    pl = pow(d, -alpha);
    gain = gamma_rv_generator(shape, scale);
}


PT vectorAddition(PT pt1, PT pt2)
{
    double x_1 = pt1.r * cos(pt1.theta);
    double y_1 = pt1.r * sin(pt1.theta);
    double x_2 = pt2.r * cos(pt2.theta);
    double y_2 = pt2.r * sin(pt2.theta);
    double x_new = x_1 + x_2;
    double y_new = y_1 + y_2;
    PT pt_new;
    pt_new.r = sqrt(pow(x_new, 2) + pow(y_new, 2));
    if(y_new >= 0) pt_new.theta = acos(x_new/pt_new.r);
    else  pt_new.theta = 2 * PI - acos(x_new/pt_new.r);
    return pt_new;
}

/**check if node from is in the main lobe of node to**/
bool isInMainLobe(Node from, Node to)
{
    double x_from = from.getPosition().r * cos(from.getPosition().theta);  ///x coordinate of node from
    double y_from = from.getPosition().r * sin(from.getPosition().theta);  ///y coordinate of node from
    double x_to = to.getPosition().r * cos(to.getPosition().theta);        ///x coordinate of node to
    double y_to = to.getPosition().r * sin(to.getPosition().theta);        ///y coordinate of node to
    double x_delta = x_from - x_to;                                        ///x coordinate of the vector that starts at node to and ends at node from
    double y_delta = y_from - y_to;                                        ///y coordinate of the above vector
    double x_boresight = cos(to.getAntenna().AOB);                         ///x coordinate of the boresight vector of node to
    double y_boresight = sin(to.getAntenna().AOB);                         ///y coordinate of the boresight vector of node to
    ///calculate the angle between the two vectors
    double delta_theta = acos((x_delta * x_boresight + y_delta * y_boresight)/sqrt(pow(x_delta, 2) + pow(y_delta, 2)));
    if (delta_theta >= 0.5 * to.getAntenna().beamwidth_m) return false;
    else return true;
}

double randomAntennaGain(Node tx, Node rx)
{
    double gain;
    /// if tx is in rx's main lobe and rx is in tx's main lobe
    if (isInMainLobe(tx, rx) && isInMainLobe(rx, tx))
        gain = tx.getAntenna().main_lobe_gain * rx.getAntenna().main_lobe_gain;
    /// if tx is in rx's main lobe and rx is in tx's back lobe
    if (isInMainLobe(tx, rx) && !isInMainLobe(rx, tx))
        gain = tx.getAntenna().back_lobe_gain * rx.getAntenna().main_lobe_gain;
    /// if tx is in rx's back lobe and rx is in tx's main lobe
    if (!isInMainLobe(tx, rx) && isInMainLobe(rx, tx))
        gain = tx.getAntenna().main_lobe_gain * rx.getAntenna().back_lobe_gain;
    /// if tx is in rx's back lobe and rx is in tx's back lobe
    if (!isInMainLobe(tx, rx) && !isInMainLobe(rx, tx))
        gain = tx.getAntenna().back_lobe_gain * rx.getAntenna().back_lobe_gain;
    return gain;
}

/*--------------------Class Node Implementation-------------------*/
Node::Node(){}

void Node::setPosition(PT ps)
{
    position = ps;
}

void Node::setPower(double tp)
{
    power = tp;
}


void Node::setDecodeTH(double dt)
{
    decodeTh = dt;
}

void Node::setSINR(double rcvP, double inter, double noise)
{
    SINR = rcvP / (inter + noise);
}

void Node::setAntenna(AT a)
{
    antenna = a;
}

double Node::transmit()
{
    return power;
}

void Node::decode()
{
    (SINR > decodeTh)?dFlag = 1:dFlag = 0;
}

PT Node::getPosition()
{
    return position;
}

AT Node::getAntenna()
{
    return antenna;
}

double Node::dis2RefNode(Node ref_node)
{
    double r1 = position.r;
    double r2 = ref_node.getPosition().r;
    double theta1 = position.theta;
    double theta2 = ref_node.getPosition().theta;

    return sqrt(pow(r1, 2) + pow(r2, 2) - 2 * r1 * r2 * cos(abs(theta1 - theta2)));
}

int Node::getDecodeFlag()
{
    return dFlag;
}

double Node::getSINR()
{
    return SINR;
}

double Node::getDecodeTH()
{
    return decodeTh;
}

Node::~Node(){}


/*------------------Class Network Implementation------------------*/
/**IO operations**/
Network::Network()
{
    transmitters = NULL;
    receivers = NULL;
    eavesdroppers = NULL;
    potentialJammers = NULL;
    index_Jammers = NULL;
    num_J = 0;
}

void Network::setRadius(double r)
{
    radius = r;
}

void Network::setDensity_D(double den_D)
{
    density_D = den_D;
    mean_D = PI * pow(radius, 2) * density_D;
}

void Network::setDensity_PJ(double den_PJ)
{
    density_PJ = den_PJ;
    mean_PJ = PI * pow(radius, 2) * density_PJ;
}

void Network::setDensity_E(double den_E)
{
    density_E = den_E;
    mean_E = PI * pow(radius, 2) * density_E;
}

void Network::set_q(double q)
{
    this->q = q;
}

int Network::getNumD()
{
    return num_D;
}

/*parameters:
tp: transmission power
w_t, w_r; main beam width of the tx and rx, resp.
g_b_t, g_b_r: back lobe gain of the tx and rx, resp.
g_m_t, g_m_r: main lobe gain of the tx and rx, resp.
dt:  decoding threshold of the rx*/
void Network::initD2DPairs(double tp, double w_t, double g_b_t, double g_m_t, double w_r, double g_b_r, double g_m_r, double dt, int flag)
{
    int i = 0;
    //deleteD2DPairs();
    /*if (0 != mean_D)
        num_D = poisson_rv_generator(mean_D);      ///generate the number of D2D pairs (poisson distributed)
    else
        num_D = 1;
    
    if(!num_D) { 
        cout<<"No nodes generated"<<endl;
        return;
    }*/
    num_D = 1 + poisson_rv_generator(mean_D);

    transmitters = new Node [num_D];
    receivers = new Node [num_D];
    PT pt_t, pt_r, temp_pt;
    for (i = 0; i < num_D; i++)
    {
        AT a_t, a_r;
        if (0 == i) 
        {       
            if(0 == flag)
            {///place the receiver at the origin
                pt_r.r = 0;
                pt_r.theta = 0;
                pt_t.r = R0;                
                pt_t.theta = uniform_rv_generator(0.0, 1.0) * 2 * PI;
                a_r.AOB = pt_t.theta;                 ///rx's antenna points towards its tx
                a_t.AOB = a_r.AOB + PI;               ///tx's antenna points towards its rx
                if(a_t.AOB > 2 * PI) a_t.AOB -= 2*PI; ///tx's AOB = (PI + rx' AOB)%(2*PI)
            }
            if(1 == flag)
            {///place the transmitter at the origin
                pt_t.r = 0;
                pt_t.theta = 0;
                pt_r.r = R0;
                pt_r.theta = uniform_rv_generator(0.0, 1.0) * 2 * PI;
                a_t.AOB = pt_r.theta;
                a_r.AOB = a_t.AOB + PI;
                if(a_r.AOB > 2 * PI) a_r.AOB -= 2*PI;
            }
        }
        else{
            pt_t.r = sqrt(uniform_rv_generator(0.0, 1.0)) * radius;
            pt_t.theta = uniform_rv_generator(0.0, 1.0) * 2 * PI;
            a_t.AOB = uniform_rv_generator(0.0, 1.0) * 2 * PI; /// tx's AOB is uniformly distributed in [0, 2*PI]
            temp_pt.r = R0;
            temp_pt.theta = a_t.AOB;              ///tx's antenna points towards its rx
            pt_r = vectorAddition(pt_t, temp_pt);
            a_r.AOB = a_t.AOB + PI;               ///rx's antenna points towards its tx
            if(a_r.AOB > 2 * PI) a_r.AOB -= 2*PI; ///rx's AOB = (PI + tx' AOB)%(2*PI)
        }
        /**initialize the i-th transmitter **/
        a_t.beamwidth_m = w_t;
        a_t.back_lobe_gain = g_b_t;
        a_t.main_lobe_gain = g_m_t;
        transmitters[i].setPosition(pt_t) ;
        transmitters[i].setPower(tp);
        transmitters[i].setAntenna(a_t);

        /**initialize the i-th receiver*/
        a_r.beamwidth_m = w_r;
        a_r.back_lobe_gain = g_b_r;
        a_r.main_lobe_gain = g_m_r;
        receivers[i].setPosition(pt_r);
        receivers[i].setAntenna(a_r);
        receivers[i].setDecodeTH(dt);
    }
}

void Network::initJammers(double tp, double w_t, double g_b_t, double g_m_t)
{
    int i = 0, j = 0;
    //deleteJammers();
    num_PJ = poisson_rv_generator(mean_PJ);      ///generate the number of D2D pairs (poisson distributed)
    potentialJammers = new Node [num_PJ];
    associatedReceivers = new int [num_PJ];
    index_Jammers = new int [num_PJ];
    num_J = 0;
    double min_d = 2 * radius;  
    PT pt_t;

    for (i = 0; i < num_PJ; i++)
    {
        /**initialize the i-th potential jammer**/
     	AT a_t;
        pt_t.r = sqrt(uniform_rv_generator(0.0, 1.0)) * radius;
        pt_t.theta = uniform_rv_generator(0.0, 1.0) * 2 * PI;
        a_t.AOB = uniform_rv_generator(0.0, 1.0) * 2 * PI; /// tx's AOB is uniformly distributed in [0, 2*PI]
        a_t.beamwidth_m = w_t;
        a_t.back_lobe_gain = g_b_t;
        a_t.main_lobe_gain = g_m_t;
        potentialJammers[i].setPosition(pt_t) ;
        potentialJammers[i].setPower(tp);
        potentialJammers[i].setAntenna(a_t);

        /**LDCJ scheme**/
        min_d = 2 * radius;
        for(j = 0; j < num_D; j++)
        {
            double d = potentialJammers[i].dis2RefNode(receivers[j]);
            if(d < min_d) {
               min_d = d;
               associatedReceivers[i] = j;
            }
        }
        
        if(min_d <= R_LOS &&  uniform_rv_generator(0.0, 1.0) < q * (1 - P_LOS) || min_d > R_LOS){
           index_Jammers[num_J++] = i;
        }
        /**End of SCJ scheme**/
    }
}


void Network::initEvesdroppers(double w_e, double g_b_e, double g_m_e, double dt)
{
    //deleteEavedroppers();
    int i = 0;
    num_E = poisson_rv_generator(mean_E);      ///generate the number of eavesdroppers (poisson distributed)
    if(num_E) eavesdroppers = new Eavesdropper [num_E];
    PT pt_e;
    for (i = 0; i < num_E; i++)
    {
        pt_e.r = sqrt(uniform_rv_generator(0.0, 1.0)) * radius;
        pt_e.theta = uniform_rv_generator(0.0, 1.0) * 2 * PI;
        AT a_e;
        a_e.AOB = uniform_rv_generator(0.0, 1.0) * 2 * PI;
        a_e.back_lobe_gain = g_b_e;
        a_e.beamwidth_m = w_e;
        a_e.main_lobe_gain = g_m_e;
        eavesdroppers[i].setPosition(pt_e) ;	  ///uniformly distribute all the eavesdroppers in the network
        eavesdroppers[i].setDecodeTH(dt);
        eavesdroppers[i].setAntenna(a_e);
    }
}


/**interfers: indices of interfers**/
double Network::interference2RefNode(int index, int type_interfer, int type_receiver)
{
    int i = 0;
    double total_inteference = 0.0;
    double channel_gain = 0.0;
    double path_loss = 0;
    double antenna_gain = 0;
    Node ref_node;
    double d = 0;

    if(RECEIVER == type_receiver) ref_node = receivers[index];
    else ref_node = eavesdroppers[index]; 
    
    if(TRANSMITTER == type_interfer)
    {
        for (i = 0; i < num_D; i++)
        {
            if(i != index){//typical pair: index 0
                d = transmitters[i].dis2RefNode(ref_node);
                if(d < R_LOS){
                   randomPathLossAndChannelGain(d, path_loss, channel_gain);
                   antenna_gain = randomAntennaGain(transmitters[i], ref_node);
                   total_inteference += transmitters[i].transmit() * channel_gain * path_loss * antenna_gain;
                }
            }
        }
    }
   
    if(JAMMER == type_interfer)
    {
       for (int i = 0; i < num_J; i++)
       {
           d = potentialJammers[index_Jammers[i]].dis2RefNode(ref_node);
           if(d < R_LOS)
           {
               if (RECEIVER == type_receiver && associatedReceivers[index_Jammers[i]] == index)
               { //the jammer is associated with the receiver with index index
                   path_loss = pow(d, -ALPHA_N);
                   channel_gain = gamma_rv_generator((double)NN, (double)1/NN);
               }
               else //the jammer is not associated with the receiver index or the receiving node is an eavesdropper
               {
                   randomPathLossAndChannelGain(d, path_loss, channel_gain);
               }
               antenna_gain = randomAntennaGain(potentialJammers[index_Jammers[i]], ref_node);
               total_inteference += potentialJammers[index_Jammers[i]].transmit() * channel_gain * path_loss * antenna_gain;
           }
       }
    }

    return total_inteference;
}

bool Network::legitmateTransmission()
{
    double rcvP, inter_T = 0.0, inter_J = 0.0;
    double transmit_power, channel_gain, antenna_gain, path_loss;
    double d = receivers[0].dis2RefNode(transmitters[0]);
    transmit_power =  transmitters[0].transmit();
    channel_gain = gamma_rv_generator((double)NL, (double)1/NL);
    antenna_gain = transmitters[0].getAntenna().main_lobe_gain * receivers[0].getAntenna().main_lobe_gain;
    path_loss = pow(d, -ALPHA_L);
    rcvP = transmit_power * channel_gain * antenna_gain * path_loss;
    inter_T = interference2RefNode(0, TRANSMITTER, RECEIVER);
    inter_J = interference2RefNode(0, JAMMER, RECEIVER);
    receivers[0].setSINR(rcvP, inter_T + inter_J, NOISE);
    receivers[0].decode();
    return receivers[0].getDecodeFlag();
}

bool Network::eavesdropping()
{
    double rcvP, inter_T, inter_J;
    double transmit_power, channel_gain, antenna_gain, path_loss;
    transmit_power =  transmitters[0].transmit();
    bool flag = false;
    for(int i = 0; i < num_E; i++)
    {
        double d = eavesdroppers[i].dis2RefNode(transmitters[0]);
        randomPathLossAndChannelGain(d, path_loss, channel_gain);
        antenna_gain = randomAntennaGain(transmitters[0], eavesdroppers[i]);
        rcvP = transmit_power * channel_gain * antenna_gain * path_loss;
        inter_T = interference2RefNode(i, TRANSMITTER, EAVESDROPPER);
        inter_J = interference2RefNode(i, JAMMER, EAVESDROPPER);
        eavesdroppers[i].setSINR(rcvP, inter_T + inter_J, NOISE);
        eavesdroppers[i].decode();
        if (eavesdroppers[i].getDecodeFlag())
        {
            flag = true;
            break;
        }
    }
    return flag;
}


void Network::deleteD2DPairs()
{
    if(transmitters) delete [] transmitters;
    if(receivers) delete [] receivers;
}


void Network::deleteJammers()
{
   if(potentialJammers) delete [] potentialJammers;
   if(index_Jammers) delete [] index_Jammers;
}


void Network::deleteEavedroppers()
{
    if(eavesdroppers) delete [] eavesdroppers;
}


/*-------------------Class Simulator Implementation------------------*/

Simulator::Simulator()
{
}

void Simulator::setLoopNum(int ln)
{
    loop_num = ln;
}

void Simulator::initParameters(Parameters pm)
{
    pmts = pm;
}

double Simulator::connectionProbability()
{
    int i = 0, cnt = 0;
    double res = 0.0;
    netW.setRadius(pmts.NW_R);
    netW.setDensity_D(pmts.den_d);
    netW.setDensity_PJ(pmts.den_pj);
    netW.setDensity_E(pmts.den_e);
    netW.set_q(pmts.q);
    int loop_actul = 0;

    for (i = 0; i < loop_num; i++)
    {
        netW.initD2DPairs(pmts.tp, pmts.w_t, pmts.g_b_t, pmts.g_m_t, pmts.w_r, pmts.g_b_r, pmts.g_m_r, pmts.dt_d, 0);
        if(netW.getNumD()){
            //loop_actul++;
            netW.initJammers(pmts.tp, pmts.w_t, pmts.g_b_t, pmts.g_m_t);
            if(netW.legitmateTransmission()) cnt++;
        }
    }
    res = (double)cnt/loop_num;
    return res;
}

double Simulator::secrecyProbabiltiy()
{
    int i = 0, cnt = 0;
    double res = 0.0;
    int loop_actul = 0;
    netW.setRadius(pmts.NW_R);
    netW.setDensity_D(pmts.den_d);
    netW.setDensity_PJ(pmts.den_pj);
    netW.setDensity_E(pmts.den_e);
    netW.set_q(pmts.q);

    for (i = 0; i < loop_num; i++)
    {
        netW.initD2DPairs(pmts.tp, pmts.w_t, pmts.g_b_t, pmts.g_m_t, pmts.w_r, pmts.g_b_r, pmts.g_m_r, pmts.dt_d, 1);
        netW.initJammers(pmts.tp, pmts.w_t, pmts.g_b_t, pmts.g_m_t);
        netW.initEvesdroppers(pmts.w_e, pmts.g_b_e, pmts.g_m_e, pmts.dt_e);
        if(!netW.eavesdropping()) cnt++;
    }
    return (double)cnt/loop_num;
}

Simulator::~Simulator()
{
}
