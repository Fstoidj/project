import javax.swing.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
public class project {


    static ArrayList<node> N=new ArrayList<node>(0);
    static ArrayList<nodes> U=new ArrayList<nodes>(0);
    static ArrayList<resistor> R=new ArrayList<resistor>(0);
    static ArrayList<currentSource> CS=new ArrayList<currentSource>(0);
    static ArrayList<voltageSource> VS=new ArrayList<voltageSource>(0);
    static ArrayList<capacitor> C=new ArrayList<capacitor>(0);
    static ArrayList<inductor> L=new ArrayList<inductor>(0);


    public static int searchNode(String s){
        for(int i=0;i<N.size();i++){
            if(N.get(i).name.equals(s)){
                return i;
            }
        }
        return -1;
    }
    public static int searchUnion(int u){
        for(int i=0;i<U.size();i++){
            if(U.get(i).union==u){
                return i;
            }
        }
        return -1;
    }


    public static class node{
        String name, output;
        double volt2=0, volt1=0;
        int union;
        node(String s){
            name=s;
        }
        public int equal(node n){
            if(n.name.equals(name)){
                return 1;
            }
            return -1;
        }
    }
    public static class nodes{
        int union;
        ArrayList<node> n=new ArrayList<node>(0);
    }
    abstract public static class branch{
        String output;
        node n1, n2;
        double I2=0, I1=0;
    }
    public static class resistor extends branch {
        String NameR;
        double R;
        resistor(String s){
            NameR=s.substring(0,s.indexOf(" "));
            s=s.substring(s.indexOf(" ")+1);
            node n=new node(s.substring(0, s.indexOf(" ")));
            int i=searchNode(s.substring(0,s.indexOf(" ")));
            nodes nodz =new nodes();
            if(i==-1){
                n.union=N.size();
                nodz.n.add(n);
                nodz.union=n.union;
                U.add(nodz);
                N.add(n);
                n1=n;
            }
            else{
                n1=N.get(i);
            }
            s=s.substring(s.indexOf(" ")+1);
            nodz =new nodes();
            n=new node(s.substring(0, s.indexOf(" ")));
            i=searchNode(s.substring(0,s.indexOf(" ")));
            if(i==-1){
                n.union=N.size();
                nodz.n.add(n);
                nodz.union=n.union;
                U.add(nodz);
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
            if(s.indexOf("m")!=-1){
                R=0.001;
                s=s.replaceAll("m","");
            }
            else if(s.indexOf("u")!=-1){
                R=0.000001;
                s=s.replaceAll("u", "");
            }
            else if(s.indexOf("n")!=-1){
                R=0.000000001;
                s=s.replaceAll("n", "");
            }
            else if(s.indexOf("p")!=-1){
                R=0.000000000001;
                s=s.replaceAll("p", "");
            }
            else {
                R=1;
            }
            R*=Double.parseDouble(s);
        }
    }
    public static class capacitor extends branch {
        String NameC;
        double C;
        ArrayList<currentSource> connectedToN2Currents=new ArrayList<currentSource>(0);
        ArrayList<voltageSource> connectedToN2Voltages=new ArrayList<voltageSource>(0);
        capacitor(String s){
            NameC=s.substring(0,s.indexOf(" "));
            s=s.substring(s.indexOf(" ")+1);
            node n=new node(s.substring(0, s.indexOf(" ")));
            nodes nodz =new nodes();
            int i=searchNode(s.substring(0,s.indexOf(" ")));
            if(i==-1){
                n.union=N.size();
                nodz.n.add(n);
                nodz.union=n.union;
                U.add(nodz);
                N.add(n);
                n1=n;
            }
            else{
                n1=N.get(i);
            }
            s=s.substring(s.indexOf(" ")+1);
            nodz =new nodes();
            n=new node(s.substring(0, s.indexOf(" ")));
            i=searchNode(s.substring(0,s.indexOf(" ")));
            if(i==-1){
                n.union=N.size();
                nodz.n.add(n);
                nodz.union=n.union;
                U.add(nodz);
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
            if(s.indexOf("m")!=-1){
                C=0.001;
                s=s.replaceAll("m","");
            }
            else if(s.indexOf("u")!=-1){
                C=0.000001;
                s=s.replaceAll("u", "");
            }
            else if(s.indexOf("n")!=-1){
                C=0.000000001;
                s=s.replaceAll("n", "");
            }
            else if(s.indexOf("p")!=-1){
                C=0.000000000001;
                s=s.replaceAll("p", "");
            }
            else {
                C=1;
            }
            C*=Double.parseDouble(s);
        }
    }
    public static class inductor extends branch {
        String NameL;
        double L;
        inductor(String s){
            NameL=s.substring(0,s.indexOf(" "));
            s=s.substring(s.indexOf(" ")+1);
            node n=new node(s.substring(0, s.indexOf(" ")));
            nodes nodz =new nodes();
            int i=searchNode(s.substring(0,s.indexOf(" ")));
            if(i==-1){
                n.union=N.size();
                nodz.n.add(n);
                nodz.union=n.union;
                U.add(nodz);
                N.add(n);
                n1=n;
            }
            else{
                n1=N.get(i);
            }
            s=s.substring(s.indexOf(" ")+1);
            n=new node(s.substring(0, s.indexOf(" ")));
            i=searchNode(s.substring(0,s.indexOf(" ")));
            nodz =new nodes();
            if(i==-1){
                n.union=N.size();
                nodz.n.add(n);
                nodz.union=n.union;
                U.add(nodz);
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
            if(s.indexOf("m")!=-1){
                L=0.001;
                s=s.replaceAll("m","");
            }
            else if(s.indexOf("u")!=-1){
                L=0.000001;
                s=s.replaceAll("u", "");
            }
            else if(s.indexOf("n")!=-1){
                L=0.000000001;
                s=s.replaceAll("n", "");
            }
            else if(s.indexOf("p")!=-1){
                L=0.000000000001;
                s=s.replaceAll("p", "");
            }
            else {
                L=1;
            }
            L*=Double.parseDouble(s);
        }
    }
    public static class voltageSource extends branch {
        String NameV;
        double V, A, w, p;
        voltageSource(String s){
            NameV=s.substring(0,s.indexOf(" "));
            s=s.substring(s.indexOf(" ")+1);
            node n=new node(s.substring(0, s.indexOf(" ")));
            nodes nodz =new nodes();
            int i=searchNode(s.substring(0,s.indexOf(" ")));
            if(i==-1){
                n.union=N.size();
                nodz.n.add(n);
                U.add(nodz);
                N.add(n);
                n1=n;
            }
            else{
                n1=N.get(i);
            }
            s=s.substring(s.indexOf(" ")+1);
            n=new node(s.substring(0, s.indexOf(" ")));
            i=searchNode(s.substring(0,s.indexOf(" ")));
            int j=searchUnion(n1.union);
            if(i==-1){
                n.union=n1.union;
                U.get(j).n.add(n);
                N.add(n);
                n2=n;
            }
            else{
                int k=searchUnion(N.get(i).union);
                for(int l=0;l<U.get(k).n.size();l++){
                    U.get(k).n.get(k).union=U.get(j).union;
                    U.get(j).n.add(U.get(k).n.get(k));
                }
                U.remove(k);
                n2=N.get(i);
            }
            s=s.substring(s.indexOf(" ")+1);
            s=s.replaceAll("k", "000");
            s=s.replaceAll("M", "000000");
            s=s.replaceAll("G", "000000000");
            if(s.indexOf("m")!=-1){
                V=0.001;
                s=s.replaceAll("m", "");
            }
            else if(s.indexOf("u")!=-1){
                V=0.000001;
                s=s.replaceAll("u", "");
            }
            else if(s.indexOf("n")!=-1){
                V=0.000000001;
                s=s.replaceAll("n", "");
            }
            else if(s.indexOf("p")!=-1){
                V=0.000000000001;
                s=s.replaceAll("p", "");
            }
            else {
                V=1;
            }
            V*=Double.parseDouble(s.substring(0, s.indexOf(" ")));
            s=s.substring(s.indexOf(" ")+1, s.length());
            A=Double.parseDouble(s.substring(0,s.indexOf(" ")));
            s=s.substring(s.indexOf(" ")+1, s.length());
            w=2*Math.PI*Double.parseDouble(s.substring(0,s.indexOf(" ")));
            s=s.substring(s.indexOf(" ")+1, s.length());
            p=Double.parseDouble(s);
            V+=A*Math.sin(p);
        }
    }
    public static class currentSource extends branch {
        String NameI;
        double A, w, p;
        currentSource(String s){
            NameI=s.substring(0,s.indexOf(" "));
            s=s.substring(s.indexOf(" ")+1);
            node n=new node(s.substring(0, s.indexOf(" ")));
            nodes nodz=new nodes();
            int i=searchNode(s.substring(0,s.indexOf(" ")));
            if(i==-1){
                n.union=N.size();
                nodz.n.add(n);
                nodz.union=n.union;
                U.add(nodz);
                N.add(n);
                n1=n;
            }
            else{
                n1=N.get(i);
            }
            s=s.substring(s.indexOf(" ")+1);
            n=new node(s.substring(0, s.indexOf(" ")));
            nodz =new nodes();
            i=searchNode(s.substring(0,s.indexOf(" ")));
            if(i==-1){
                n.union=N.size();
                nodz.n.add(n);
                nodz.union=n.union;
                U.add(nodz);
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
            if(s.indexOf("m")!=-1){
                I1=0.001;
                s=s.replaceAll("m", "");
            }
            else if(s.indexOf("u")!=-1){
                I1=0.000001;
                s=s.replaceAll("u", "");
            }
            else if(s.indexOf("n")!=-1){
                I1=0.000000001;
                s=s.replaceAll("n", "");
            }
            else if(s.indexOf("p")!=-1){
                I1=0.000000000001;
                s=s.replaceAll("p", "");
            }
            else {
                I1=1;
            }
            I1*=Double.parseDouble(s.substring(0, s.indexOf(" ")));
            s=s.substring(s.indexOf(" ")+1, s.length());
            A=Double.parseDouble(s.substring(0,s.indexOf(" ")));
            s=s.substring(s.indexOf(" ")+1, s.length());
            w=2*Math.PI*Double.parseDouble(s.substring(0,s.indexOf(" ")));
            s=s.substring(s.indexOf(" ")+1, s.length());
            p=Double.parseDouble(s);
            I1+=A*Math.sin(p);
            I2=I1;
        }
    }
    public class diode extends branch {
        String NameD;
        void checdiod(){
            if(I1<0){
                I1=0;
            }
            else if(n1.volt1-n2.volt1>0){
                n1.volt1=n2.volt1;
            }
        }
    }




    public static void calcNodeVolts(double dT, double dV, double dI){
        double i1=0, i2=0;

        for(int i=1;i<U.size();i++){
            for(int j=0;j<U.get(i).n.size();j++){
                for(int k=0; k<R.size();k++){
                    if(R.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i1+=(R.get(k).n2.volt1-R.get(k).n1.volt1)/R.get(k).R;
                    }
                    else if(R.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i1-=(R.get(k).n2.volt1-R.get(k).n1.volt1)/R.get(k).R;
                    }
                }
                for(int k=0; k<CS.size();k++){
                    if(CS.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i1+=CS.get(k).I1;
                    }
                    else if(CS.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i1-=CS.get(k).I1;
                    }
                }
            }

            for(int j=0;j<U.get(i).n.size();j++){
                for(int k=0; k<R.size();k++){
                    if(R.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i2+=(R.get(k).n2.volt1-(R.get(k).n1.volt1+dV))/R.get(k).R;
                    }
                    else if(R.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i2-=((R.get(k).n2.volt1+dV)-R.get(k).n1.volt1)/R.get(k).R;
                    }
                }
                for(int k=0; k<CS.size();k++){
                    if(CS.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i2+=CS.get(k).I1;
                    }
                    else if(CS.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i2-=CS.get(k).I1;
                    }
                }
            }
            U.get(i).n.get(0).volt2+=((Math.abs(i1)-Math.abs(i2))*dV)/dI;
            i1=0;
            i2=0;
        }


        for(int i=0;i<U.size();i++){
            U.get(i).n.get(0).volt1=U.get(i).n.get(0).volt2;
            System.out.println("Node: "+U.get(i).n.get(0).name+" Volt is: "+U.get(i).n.get(0).volt1);
        }
    }



    public static class graphProject extends JFrame{


    }

    public static void main(String[] args) {
        node n;
        resistor r;
        currentSource cs;
        voltageSource vs;
        capacitor c;
        inductor l;
        double T = -1, dT = -1, dV=-1, dI=-1;
        Scanner sc = new Scanner(System.in);
        String s = sc.nextLine();
        s = s.trim();
        s = s.replaceAll("( )+", " ");

        while (!s.equals(".end")) {
            if (s.charAt(0) == 'R') {
                r = new resistor(s);
                R.add(r);
            }
            else if (s.charAt(0) == 'L') {
                l = new inductor(s);
                L.add(l);
            }
            else if (s.charAt(0) == 'I') {
                cs = new currentSource(s);
                CS.add(cs);
            }
            else if (s.charAt(0) == 'V') {
                vs = new voltageSource(s);
                VS.add(vs);
            }
            else if (s.charAt(0) == 'C') {
                c = new capacitor(s);
                C.add(c);
            }
            else if (s.charAt(0) == 'F') { }
            else if (s.indexOf(".tran") != -1) {
                s = s.substring(s.indexOf(" ") + 1);
                s = s.replaceAll("k", "000");
                s = s.replaceAll("M", "000000");
                s = s.replaceAll("G", "000000000");
                if (s.indexOf("m") != -1) {
                    T = 0.001;
                    s = s.replaceAll("m", "");
                } else if (s.indexOf("u") != -1) {
                    T = 0.000001;
                    s = s.replaceAll("u", "");
                } else if (s.indexOf("n") != -1) {
                    T = 0.000000001;
                    s = s.replaceAll("n", "");
                } else if (s.indexOf("p") != -1) {
                    T = 0.000000000001;
                    s = s.replaceAll("p", "");
                } else {
                    T = 1;
                }
                T *= Double.parseDouble(s);
            }
            else if (s.indexOf("dT") != -1) {
                s = s.substring(s.indexOf(" ") + 1);
                s = s.replaceAll("k", "000");
                s = s.replaceAll("M", "000000");
                s = s.replaceAll("G", "000000000");
                if (s.indexOf("m") != -1) {
                    dT = 0.001;
                    s = s.replaceAll("m", "");
                } else if (s.indexOf("u") != -1) {
                    dT = 0.000001;
                    s = s.replaceAll("u", "");
                } else if (s.indexOf("n") != -1) {
                    dT = 0.000000001;
                    s = s.replaceAll("n", "");
                } else if (s.indexOf("p") != -1) {
                    dT = 0.000000000001;
                    s = s.replaceAll("p", "");
                } else {
                    dT = 1;
                }
                dT *= Double.parseDouble(s);

            }
            else if (s.indexOf("dV") != -1) {
                s = s.substring(s.indexOf(" ") + 1);
                s = s.replaceAll("k", "000");
                s = s.replaceAll("M", "000000");
                s = s.replaceAll("G", "000000000");
                if (s.indexOf("m") != -1) {
                    dV = 0.001;
                    s = s.replaceAll("m", "");
                } else if (s.indexOf("u") != -1) {
                    dV = 0.000001;
                    s = s.replaceAll("u", "");
                } else if (s.indexOf("n") != -1) {
                    dV = 0.000000001;
                    s = s.replaceAll("n", "");
                } else if (s.indexOf("p") != -1) {
                    dV = 0.000000000001;
                    s = s.replaceAll("p", "");
                } else {
                    dV = 1;
                }
                dV *= Double.parseDouble(s);

            }
            else if (s.indexOf("dI") != -1) {
                s = s.substring(s.indexOf(" ") + 1);
                s = s.replaceAll("k", "000");
                s = s.replaceAll("M", "000000");
                s = s.replaceAll("G", "000000000");
                if (s.indexOf("m") != -1) {
                    dI = 0.001;
                    s = s.replaceAll("m", "");
                } else if (s.indexOf("u") != -1) {
                    dI = 0.000001;
                    s = s.replaceAll("u", "");
                } else if (s.indexOf("n") != -1) {
                    dI = 0.000000001;
                    s = s.replaceAll("n", "");
                } else if (s.indexOf("p") != -1) {
                    dI = 0.000000000001;
                    s = s.replaceAll("p", "");
                } else {
                    dI = 1;
                }
                dI *= Double.parseDouble(s);

            }

            s = sc.nextLine();
            s = s.trim();
            s = s.replaceAll("( )+", " ");

        }

        //testing union merging in VS
        System.out.println("\n");
        for (int j = 0; j < U.size(); j++) {
            System.out.println(U.get(j).n.get(0).name+" : "+U.get(j).n.get(0).union+" in: "+U.get(j).union);
        }


        if(dT>0&&T>0&&dI>0&&dV>0){
            for(int i=0;i<T/dT;i++) {
                System.out.println("iteration number: "+i);
                calcNodeVolts(dT, dV, dI);
                System.out.println("\n");
            }
        }


    }
}
