public class project {
    public class node{
        String name;
        double volt;
    }
    abstract public class branch{
        node n1, n2;
        double I;
    }
    public class resistor extends branch {
        String NameR;
        double R;
        void calc_I(){
            I=(n1.volt-n2.volt)/R;
        }
        void calc_v1(){
            n1.volt=n2.volt+I*R;
        }
        void calc_v2(){
            n2.volt=n1.volt-I*R;
        }
    }
    public class capacitor extends branch {
        String NameC;
        double C, Q=0, Q0=0;
        void calc_QwithV(double dt){
            Q0=Q;
            Q=C*(n1.volt-n2.volt);
            I=(Q-Q0)/dt;
        }
        void calc_Q1withI(double dt){
            Q0=Q;
            Q=Q0+dt*I;
            n1.volt=(Q-Q0)/C+n2.volt;
        }
        void calc_Q2withI(double dt){
            Q0=Q;
            Q=Q0+dt*I;
            n2.volt=n1.volt-(Q-Q0)/C;
        }
    }
    public class inductor extends branch {
        String NameL;
        double L, I0=0;
        void calcI(double dt){
            I0=I;
            I=I0+dt*(n1.volt-n2.volt)/L;
        }
        void calcV1(double dt){
            n1.volt=L*(I-I0)/dt+n2.volt;
        }
        void calcV2(double dt){
            n2.volt=n1.volt-L*(I-I0)/dt;
        }
    }
    public class voltageSource extends branch {
        String NameV;
        double V, freq;
        void calcV1(){
            n1.volt=V+n2.volt;
        }
        void calcV2(){
            n2.volt=n1.volt-V;
        }
    }
    public class currentSource extends branch {
        String NameI;
        double freq;
    }
    public class diode extends branch {
        String NameD;
        void checdiod(){
            if(I<0){
                I=0;
            }
            else if(n1.volt-n2.volt>0){
                n1.volt=n2.volt;
            }
        }
    }
    public static void main(String[] args){
        voltageSource v;

    }
}
