import javax.swing.*;
import java.awt.*;
import java.awt.image.ImageObserver;
import java.awt.image.ImageProducer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
public class project {


    static ArrayList<String> input=new ArrayList<String>(0);


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
        String name;
        ArrayList<Double> outputVolt=new ArrayList<Double>(1);
        double volt=0;
        int union;
        node(String s){
            ArrayList<Double> outputVolt=new ArrayList<Double>(1);
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
        ArrayList<Double> outputCurrent=new ArrayList<Double>(0);
        node n1, n2;
        double I=0;
        public void output(){
            for(int i=0;i<n1.outputVolt.size();i++){
                System.out.print(" "+(n1.outputVolt.get(i)-n2.outputVolt.get(i))+" "+outputCurrent.get(i)+" "+(outputCurrent.get(i)*(n1.outputVolt.get(i)-n2.outputVolt.get(i))));
            }
        }
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
                I=0.001;
                s=s.replaceAll("m", "");
            }
            else if(s.indexOf("u")!=-1){
                I=0.000001;
                s=s.replaceAll("u", "");
            }
            else if(s.indexOf("n")!=-1){
                I=0.000000001;
                s=s.replaceAll("n", "");
            }
            else if(s.indexOf("p")!=-1){
                I=0.000000000001;
                s=s.replaceAll("p", "");
            }
            else {
                I=1;
            }
            I*=Double.parseDouble(s.substring(0, s.indexOf(" ")));
            s=s.substring(s.indexOf(" ")+1, s.length());
            A=Double.parseDouble(s.substring(0,s.indexOf(" ")));
            s=s.substring(s.indexOf(" ")+1, s.length());
            w=2*Math.PI*Double.parseDouble(s.substring(0,s.indexOf(" ")));
            s=s.substring(s.indexOf(" ")+1, s.length());
            p=Double.parseDouble(s);
            I+=A*Math.sin(p);
            outputCurrent.add(I);
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




    public static void calcNodeVolts(double dT, double dV, double dI){
        double i1=0, i2=0;

        U.get(0).n.get(0).outputVolt.add(0.0);
        for(int i=1;i<U.size();i++){
            for(int j=0;j<U.get(i).n.size();j++){
                for(int k=0; k<R.size();k++){
                    if(R.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i1+=(R.get(k).n2.volt-R.get(k).n1.volt)/R.get(k).R;
                    }
                    else if(R.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i1-=(R.get(k).n2.volt-R.get(k).n1.volt)/R.get(k).R;
                    }
                }
                for(int k=0; k<CS.size();k++){
                    if(CS.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i1+=CS.get(k).I;
                    }
                    else if(CS.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i1-=CS.get(k).I;
                    }
                }
            }

            for(int j=0;j<U.get(i).n.size();j++){
                for(int k=0; k<R.size();k++){
                    if(R.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i2+=(R.get(k).n2.volt-(R.get(k).n1.volt+dV))/R.get(k).R;
                    }
                    else if(R.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i2-=((R.get(k).n2.volt+dV)-R.get(k).n1.volt)/R.get(k).R;
                    }
                }
                for(int k=0; k<CS.size();k++){
                    if(CS.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i2+=CS.get(k).I;
                    }
                    else if(CS.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i2-=CS.get(k).I;
                    }
                }
            }
            U.get(i).n.get(0).outputVolt.add(U.get(i).n.get(0).volt+((Math.abs(i1)-Math.abs(i2))*dV)/dI);
            i1=0;
            i2=0;
        }

        for(int i=0, j;i<U.size();i++){
            j=U.get(i).n.get(0).outputVolt.size()-1;
            U.get(i).n.get(0).volt=U.get(i).n.get(0).outputVolt.get(j);
        }
    }
    public static void calcBranchCurrents(double dT, double dV, double dI){
        for(int i=0;i<R.size();i++){
            R.get(i).outputCurrent.add((R.get(i).n1.volt-R.get(i).n2.volt)/R.get(i).R);
        }
        for (int i=0;i<CS.size();i++){
            CS.get(i).outputCurrent.add(CS.get(i).I);
        }
    }

    public static void chapOutput(){
        for(int i=0;i<N.size();i++){
            System.out.print(U.get(i).n.get(0).name+" :");
            for (int j=0;j<U.get(0).n.get(0).outputVolt.size();j++){
                System.out.print(" "+U.get(i).n.get(0).outputVolt.get(j));
            }
            System.out.println();
        }

        for(int i=0, iR=0, iI=0, iC=0, iL=0, iV=0; i<input.size();i++){
            if(input.get(i).equals("R")){
                System.out.print(R.get(iR).NameR+" :");
                R.get(iR).output();
                iR++;
            }
            else if(input.get(i).equals("I")){
                System.out.print(CS.get(iI).NameI+" :");
                CS.get(iI).output();
                iI++;
            }
            else if(input.get(i).equals("C")){
                System.out.print(C.get(iC).NameC+" :");
                C.get(iC).output();
                iC++;
            }
            else if(input.get(i).equals("L")){
                System.out.print(L.get(iL).NameL+" :");
                L.get(iL).output();
                iL++;
            }
            else if(input.get(i).equals("V")){
                System.out.print(VS.get(iV).NameV+" :");
                VS.get(iV).output();
                iV++;
            }
            System.out.println();
        }
    }


    public static class graphProject extends JFrame{

        graphProject() {
            JLabel label=new JLabel(new ImageIcon("7878.png"));
            label.setBounds(500,500,180,140);
            add(label);

            setTitle("Circuit Graph");
            setSize(2000, 1000);
            setLayout(null);
            setVisible(true);
        }

    }

    public static void main(String[] args) {
        graphProject GraphProject=new graphProject();

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
                r.outputCurrent.add(0.0);
                R.add(r);
                input.add("R");
            }
            else if (s.charAt(0) == 'L') {
                l = new inductor(s);
                L.add(l);
                input.add("L");
            }
            else if (s.charAt(0) == 'I') {
                cs = new currentSource(s);
                CS.add(cs);
                input.add("I");
            }
            else if (s.charAt(0) == 'V') {
                vs = new voltageSource(s);
                VS.add(vs);
                input.add("V");
            }
            else if (s.charAt(0) == 'C') {
                c = new capacitor(s);
                C.add(c);
                input.add("C");
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
        for (int j = 0; j < U.size(); j++) {
            U.get(j).n.get(0).outputVolt.add(0.0);
            //System.out.println(U.get(j).n.get(0).name+" : "+U.get(j).n.get(0).union+" in: "+U.get(j).union);
        }

        if(dT>0&&T>0&&dI>0&&dV>0){
            for(int i=0;i<T/dT;i+=1) {
                calcNodeVolts(dT, dV, dI);
                calcBranchCurrents(dT, dV, dI);
            }
        }

        chapOutput();

    }
}
