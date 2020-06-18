import java.util.ArrayList;
import java.util.Scanner;

public class project {

    static ArrayList<node> N=new ArrayList<node>(0);
    static ArrayList<resistor> R=new ArrayList<resistor>(0);
    static ArrayList<currentSource> CS=new ArrayList<currentSource>(0);
    static ArrayList<voltageSource> VS=new ArrayList<voltageSource>(0);

    public static int searchNode(String s){
        for(int i=0;i<N.size();i++){
            if(N.get(i).name.equals(s)){
                return i;
            }
        }
        return -1;
    }


    public static class node{
        String name;
        double volt;
        node(String s){
            name=s;
        }
    }
    abstract public static class branch{
        node n1, n2;
        double I;
    }
    public static class resistor extends branch {
        String NameR;
        double R;
        resistor(String s){
            NameR=s.substring(0,s.indexOf(" "));
            s=s.substring(s.indexOf(" ")+1);
            node n=new node(s.substring(0, s.indexOf(" ")));
            int i=searchNode(s.substring(0,s.indexOf(" ")));
            if(i==-1){
                N.add(n);
                n1=n;
            }
            else{
                n1=N.get(i);
            }
            s=s.substring(s.indexOf(" ")+1);
            n=new node(s.substring(0, s.indexOf(" ")));
            i=searchNode(s.substring(0,s.indexOf(" ")));
            if(i==-1){
                N.add(n);
                n2=n;
            }
            else{
                n2=N.get(i);
            }
            s=s.substring(s.indexOf(" ")+1);
            s=s.replaceAll("k", "000");
            s=s.replaceAll("M", "000000");
            s=s.replaceAll("G", "000000000");
            if(s.indexOf("u")!=-1){
                R=0.001;
                s=s.replaceAll("u", "");
            }
            else if(s.indexOf("n")!=-1){
                R=0.000001;
                s=s.replaceAll("u", "");
            }
            else if(s.indexOf("p")!=-1){
                R=0.000000001;
                s=s.replaceAll("u", "");
            }
            else {
                R=1;
            }
            R*=Double.parseDouble(s);
        }
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
    public static class voltageSource extends branch {
        String NameV;
        double V, freq;
        voltageSource(String s){
            NameV=s.substring(0,s.indexOf(" "));
            s=s.substring(s.indexOf(" ")+1);
            node n=new node(s.substring(0, s.indexOf(" ")));
            int i=searchNode(s.substring(0,s.indexOf(" ")));
            if(i==-1){
                N.add(n);
                n1=n;
            }
            else{
                n1=N.get(i);
            }
            s=s.substring(s.indexOf(" ")+1);
            n=new node(s.substring(0, s.indexOf(" ")));
            i=searchNode(s.substring(0,s.indexOf(" ")));
            if(i==-1){
                N.add(n);
                n2=n;
            }
            else{
                n2=N.get(i);
            }
            s=s.substring(s.indexOf(" ")+1);
            s=s.replaceAll("k", "000");
            s=s.replaceAll("M", "000000");
            s=s.replaceAll("G", "000000000");
            if(s.indexOf("u")!=-1){
                V=0.001;
                s=s.replaceAll("u", "");
            }
            else if(s.indexOf("n")!=-1){
                V=0.000001;
                s=s.replaceAll("u", "");
            }
            else if(s.indexOf("p")!=-1){
                V=0.000000001;
                s=s.replaceAll("u", "");
            }
            else {
                V=1;
            }
            V*=Double.parseDouble(s.substring(0, s.indexOf(" ")));
        }
        void calcV1(){
            n1.volt=V+n2.volt;
        }
        void calcV2(){
            n2.volt=n1.volt-V;
        }
    }
    public static class currentSource extends branch {
        String NameI;
        double freq;
        currentSource(String s){
            NameI=s.substring(0,s.indexOf(" "));
            s=s.substring(s.indexOf(" ")+1);
            node n=new node(s.substring(0, s.indexOf(" ")));
            int i=searchNode(s.substring(0,s.indexOf(" ")));
            if(i==-1){
                N.add(n);
                n1=n;
            }
            else{
                n1=N.get(i);
            }
            s=s.substring(s.indexOf(" ")+1);
            n=new node(s.substring(0, s.indexOf(" ")));
            i=searchNode(s.substring(0,s.indexOf(" ")));
            if(i==-1){
                N.add(n);
                n2=n;
            }
            else{
                n2=N.get(i);
            }
            s=s.substring(s.indexOf(" ")+1);
            s=s.replaceAll("k", "000");
            s=s.replaceAll("M", "000000");
            s=s.replaceAll("G", "000000000");
            if(s.indexOf("u")!=-1){
                I=0.001;
                s=s.replaceAll("u", "");
            }
            else if(s.indexOf("n")!=-1){
                I=0.000001;
                s=s.replaceAll("u", "");
            }
            else if(s.indexOf("p")!=-1){
                I=0.000000001;
                s=s.replaceAll("u", "");
            }
            else {
                I=1;
            }
            I*=Double.parseDouble(s.substring(0, s.indexOf(" ")));
        }
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
        resistor r;
        currentSource cs;
        voltageSource vs;
        Scanner sc=new Scanner(System.in);
        String s=sc.nextLine();
        s=s.trim();
        s=s.replaceAll("( )+", " ");

        while (!s.equals(".end")){
            if(s.charAt(0)=='R'){
                r=new resistor(s);
                R.add(r);
            }
            else if(s.charAt(0)=='I'){
                cs=new currentSource(s);
                CS.add(cs);
            }
            else if(s.charAt(0)=='V'){
                vs=new voltageSource(s);
                VS.add(vs);
            }
            s=sc.nextLine();
            s=s.trim();
            s=s.replaceAll("( )+", " ");
        }
    }
}
