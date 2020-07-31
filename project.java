
import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.geom.Line2D;
import java.awt.image.BufferedImage;
import java.awt.image.ImageObserver;
import java.awt.image.ImageProducer;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
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
    static ArrayList<VCCS> G=new ArrayList<VCCS>(0);
    static ArrayList<CCCS> F=new ArrayList<CCCS>(0);


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


    public static ArrayList<Double> findBranchCurrent(String s){
        ArrayList<Double> a=new ArrayList<Double>(0);
        if (s.charAt(0) == 'R') {
            for (int i=0;i<R.size();i++){
                if(R.get(i).NameR.equals(s)){
                    return R.get(i).outputCurrent;
                }
            }
        }
        else if (s.charAt(0) == 'L') {
            for (int i=0;i<L.size();i++){
                if(L.get(i).NameL.equals(s)){
                    return L.get(i).outputCurrent;
                }
            }
        }
        else if (s.charAt(0) == 'I') {
            for (int i=0;i<CS.size();i++){
                if(CS.get(i).NameI.equals(s)){
                    return CS.get(i).outputCurrent;
                }
            }
        }
        else if (s.charAt(0) == 'V') {
            for (int i=0;i<VS.size();i++){
                if(VS.get(i).NameV.equals(s)){
                    return VS.get(i).outputCurrent;
                }
            }
        }
        else if (s.charAt(0) == 'C') {
            for (int i=0;i<C.size();i++){
                if(C.get(i).NameC.equals(s)){
                    return C.get(i).outputCurrent;
                }
            }
        }
        else if (s.charAt(0) == 'G') {
            for (int i=0;i<G.size();i++){
                if(G.get(i).NameG.equals(s)){
                    return G.get(i).outputCurrent;
                }
            }
        }
        else if (s.charAt(0) == 'F') {
            for (int i=0;i<F.size();i++){
                if(F.get(i).NameF.equals(s)){
                    return F.get(i).outputCurrent;
                }
            }
        }
        return a;
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
        public double moshtaghVolt(double dT, int iteration){
            double VDot=0;
            /*if(outputVolt.size()>iteration+2) {
                System.out.println(iteration);
                VDot = (outputVolt.get(iteration+1) - outputVolt.get(iteration)) / dT;
                return VDot;
            }*/
            if(iteration>0){
                VDot = (outputVolt.get(iteration) - outputVolt.get(iteration-1)) / dT;
                return VDot;
            }
            return 0;
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
            if(s.substring(0, s.indexOf(" ")).indexOf("m")!=-1){
                I=0.001;
                s=s.replaceAll("m", "");
            }
            else if(s.substring(0, s.indexOf(" ")).indexOf("u")!=-1){
                I=0.000001;
                s=s.replaceAll("u", "");
            }
            else if(s.substring(0, s.indexOf(" ")).indexOf("n")!=-1){
                I=0.000000001;
                s=s.replaceAll("n", "");
            }
            else if(s.substring(0, s.indexOf(" ")).indexOf("p")!=-1){
                I=0.000000000001;
                s=s.replaceAll("p", "");
            }
            else {
                I=1;
            }
            I*=Double.parseDouble(s.substring(0, s.indexOf(" ")));
            s=s.substring(s.indexOf(" ")+1, s.length());
            if(s.substring(0, s.indexOf(" ")).indexOf("m")!=-1){
                A=0.001;
                s=s.replaceAll("m", "");
            }
            else if(s.substring(0, s.indexOf(" ")).indexOf("u")!=-1){
                A=0.000001;
                s=s.replaceAll("u", "");
            }
            else if(s.substring(0, s.indexOf(" ")).indexOf("n")!=-1){
                A=0.000000001;
                s=s.replaceAll("n", "");
            }
            else if(s.substring(0, s.indexOf(" ")).indexOf("p")!=-1){
                A=0.000000000001;
                s=s.replaceAll("p", "");
            }
            else {
                A=1;
            }
            A*=Double.parseDouble(s.substring(0,s.indexOf(" ")));
            s=s.substring(s.indexOf(" ")+1, s.length());
            if(s.substring(0, s.indexOf(" ")).indexOf("m")!=-1){
                w=0.001;
                s=s.replaceAll("m", "");
            }
            else if(s.substring(0, s.indexOf(" ")).indexOf("u")!=-1){
                w=0.000001;
                s=s.replaceAll("u", "");
            }
            else if(s.substring(0, s.indexOf(" ")).indexOf("n")!=-1){
                w=0.000000001;
                s=s.replaceAll("n", "");
            }
            else if(s.substring(0, s.indexOf(" ")).indexOf("p")!=-1){
                w=0.000000000001;
                s=s.replaceAll("p", "");
            }
            else {
                w=1;
            }
            w*=2*Math.PI*Double.parseDouble(s.substring(0,s.indexOf(" ")));
            s=s.substring(s.indexOf(" ")+1, s.length());
            if(s.indexOf("m")!=-1){
                p=0.001;
                s=s.replaceAll("m", "");
            }
            else if(s.indexOf("u")!=-1){
                p=0.000001;
                s=s.replaceAll("u", "");
            }
            else if(s.indexOf("n")!=-1){
                p=0.000000001;
                s=s.replaceAll("n", "");
            }
            else if(s.indexOf("p")!=-1){
                p=0.000000000001;
                s=s.replaceAll("p", "");
            }
            else {
                p=1;
            }
            p*=Double.parseDouble(s);
            outputCurrent.add(I+A*Math.sin(p));
        }
    }
    public static class VCCS extends branch {
        String NameG;
        node n3, n4;
        double a;
        VCCS(String s){
            NameG=s.substring(0,s.indexOf(" "));
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
            n=new node(s.substring(0, s.indexOf(" ")));
            i=searchNode(s.substring(0,s.indexOf(" ")));
            if(i==-1){
            }
            else{
                n3=N.get(i);
            }
            s=s.substring(s.indexOf(" ")+1);
            n=new node(s.substring(0, s.indexOf(" ")));
            i=searchNode(s.substring(0,s.indexOf(" ")));
            if(i==-1){
            }
            else{
                n4=N.get(i);
            }
            s=s.replaceAll("k", "000");
            s=s.replaceAll("M", "000000");
            s=s.replaceAll("G", "000000000");
            if(s.indexOf("m")!=-1){
                a=0.001;
                s=s.replaceAll("m", "");
            }
            else if(s.indexOf("u")!=-1){
                a=0.000001;
                s=s.replaceAll("u", "");
            }
            else if(s.indexOf("n")!=-1){
                a=0.000000001;
                s=s.replaceAll("n", "");
            }
            else if(s.indexOf("p")!=-1){
                a=0.000000000001;
                s=s.replaceAll("p", "");
            }
            else {
                a=1;
            }
            a*=Double.parseDouble(s);
            I=n3.outputVolt.get(0)-n4.outputVolt.get(0);
            I*=a;
            outputCurrent.add(I);
        }
    }
    public static class CCCS extends branch {
        String NameF;
        String type;
        int index;
        ArrayList<Double> inputCurrent=new ArrayList<Double>(0);
        double a;
        CCCS(String s){
            NameF=s.substring(0,s.indexOf(" "));
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
            inputCurrent=findBranchCurrent(s.substring(0, s.indexOf(" ")));
            s=s.substring(s.indexOf(" ")+1);
            s=s.replaceAll("k", "000");
            s=s.replaceAll("M", "000000");
            s=s.replaceAll("G", "000000000");
            if(s.indexOf("m")!=-1){
                a=0.001;
                s=s.replaceAll("m", "");
            }
            else if(s.indexOf("u")!=-1){
                a=0.000001;
                s=s.replaceAll("u", "");
            }
            else if(s.indexOf("n")!=-1){
                a=0.000000001;
                s=s.replaceAll("n", "");
            }
            else if(s.indexOf("p")!=-1){
                a=0.000000000001;
                s=s.replaceAll("p", "");
            }
            else {
                a=1;
            }
            a*=Double.parseDouble(s);
            I=inputCurrent.get(0);
            I*=a;
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




    public static void calcNodeVolts(double dT, double dV, double dI, int iteraroin){
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
                        i1+=CS.get(k).outputCurrent.get(iteraroin+1);
                    }
                    else if(CS.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i1-=CS.get(k).outputCurrent.get(iteraroin+1);
                    }
                }
                for(int k=0; k<C.size();k++){
                    if(C.get(k).n1.name.equals(U.get(i).n.get(j).name)) {
                        i1 -= C.get(k).C * (C.get(k).n1.moshtaghVolt(dT,iteraroin) - C.get(k).n2.moshtaghVolt(dT,iteraroin));
                    }
                    else if(C.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i1 += C.get(k).C * (C.get(k).n1.moshtaghVolt(dT, iteraroin) - C.get(k).n2.moshtaghVolt(dT, iteraroin));
                    }
                }
                for(int k=0; k<L.size();k++){
                    if(L.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i1+=L.get(k).outputCurrent.get(iteraroin+1);
                    }
                    else if(L.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i1-=L.get(k).outputCurrent.get(iteraroin+1);
                    }
                }
                for(int k=0; k<G.size();k++){
                    if(G.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i1+=G.get(k).outputCurrent.get(iteraroin+1);
                    }
                    else if(G.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i1-=G.get(k).outputCurrent.get(iteraroin+ 1);
                    }
                }
                for(int k=0; k<F.size();k++){
                    if(F.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i1+=F.get(k).outputCurrent.get(iteraroin+ 1);
                    }
                    else if(F.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i1-=F.get(k).outputCurrent.get(iteraroin+ 1);
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
                        i2+=CS.get(k).outputCurrent.get(iteraroin+1);
                    }
                    else if(CS.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i2-=CS.get(k).outputCurrent.get(iteraroin+1);
                    }
                }
                for(int k=0; k<C.size();k++){
                    if(C.get(k).n1.name.equals(U.get(i).n.get(j).name)) {
                        i2 -= C.get(k).C * (C.get(k).n1.moshtaghVolt(dT, iteraroin)+dV/dT - C.get(k).n2.moshtaghVolt(dT, iteraroin));
                    }
                    else if(C.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i2 += C.get(k).C * (C.get(k).n1.moshtaghVolt(dT, iteraroin) - C.get(k).n2.moshtaghVolt(dT, iteraroin)-dV/dT);
                    }
                }
                for(int k=0; k<L.size();k++){
                    if(L.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i2+=L.get(k).outputCurrent.get(iteraroin+1)+dV*dT/L.get(k).L;
                    }
                    else if(L.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i2-=L.get(k).outputCurrent.get(iteraroin+ 1)-dV*dT/L.get(k).L;
                    }
                }
                for(int k=0; k<G.size();k++){
                    if(G.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i2+=G.get(k).outputCurrent.get(iteraroin+1);
                        if(G.get(k).n3.name.equals(U.get(i).n.get(j).name)){
                            i2+=G.get(k).a*dV;
                        }
                        else if(G.get(k).n4.name.equals(U.get(i).n.get(j).name)){
                            i2-=G.get(k).a*dV;
                        }
                    }
                    else if(G.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i2-=G.get(k).outputCurrent.get(iteraroin+ 1);
                        if(G.get(k).n3.name.equals(U.get(i).n.get(j).name)){
                            i2+=G.get(k).a*dV;
                        }
                        else if(G.get(k).n4.name.equals(U.get(i).n.get(j).name)){
                            i2-=G.get(k).a*dV;
                        }
                    }
                }
                for(int k=0; k<F.size();k++){
                    if(F.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i2+=F.get(k).outputCurrent.get(iteraroin+1);
                    }
                    else if(F.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i2-=F.get(k).outputCurrent.get(iteraroin+ 1);
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
    public static void calcBranchCurrents(double dT, double dV, double dI, int iteration){
        double I=0;
        for(int i=0;i<R.size();i++){
            R.get(i).outputCurrent.add((R.get(i).n1.volt-R.get(i).n2.volt)/R.get(i).R);
        }
        for (int i=0;i<CS.size();i++){
            CS.get(i).outputCurrent.add(CS.get(i).I+CS.get(i).A*Math.sin(CS.get(i).p+CS.get(i).w*(iteration+1)*dT));
        }
        for (int i=0; i<C.size();i++){
            C.get(i).outputCurrent.add(C.get(i).C * (C.get(i).n1.moshtaghVolt(dT, iteration) - C.get(i).n2.moshtaghVolt(dT, iteration)));
        }
        for (int i=0; i<L.size();i++){
            I=L.get(i).outputCurrent.get(L.get(i).outputCurrent.size()-1);
            L.get(i).outputCurrent.add(I+(dT*(L.get(i).n2.outputVolt.get(L.get(i).n2.outputVolt.size()-1)-L.get(i).n1.outputVolt.get(L.get(i).n1.outputVolt.size()-1)))/L.get(i).L);
        }
        for (int i=0;i<G.size();i++){
            G.get(i).outputCurrent.add((G.get(i).n3.outputVolt.get(0)-G.get(i).n4.outputVolt.get(0))*G.get(i).a);
        }
        for (int i=0;i<F.size();i++){
            F.get(i).outputCurrent.add(F.get(i).inputCurrent.get(iteration+1));
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

        for(int i=0, iR=0, iI=0, iC=0, iL=0, iV=0, iG=0, iF=0; i<input.size();i++){
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
            else if(input.get(i).equals("G")){
                System.out.print(G.get(iG).NameG+" :");
                G.get(iG).output();
                iG++;
            }
            else if(input.get(i).equals("F")){
                System.out.print(F.get(iF).NameF+" :");
                F.get(iF).output();
                iF++;
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

        ArrayList<Line2D> lines = new ArrayList<>();
        HashMap<String, Integer> paralel = new HashMap<>();

        graphProject() {
            //JLabel label=new JLabel(new ImageIcon("7878.png"));
            //label.setBounds(200,350,180,100);
            //add(label);
            setTitle("Circuit Graph");
            setSize(2000, 1000);
            setLayout(null);
            setVisible(true);
        }

        @Override
        public void paint(Graphics g) {
            super.paint(g);
            g.setColor(Color.BLACK);
            g.drawLine(0,0,30,1000);
            for (resistor resistor1: R) {
                drawResistor(resistor1);
            }
            for (currentSource currentSource1: CS){
                drawCS(currentSource1);
            }
            for (voltageSource voltageSource1: VS){
                drawVS(voltageSource1);
            }
            for (capacitor capacitor1: C){
                drawC(capacitor1);
            }
            for (inductor inductor1: L){
                drawL(inductor1);
            }
            drawLines();
        }

        private void drawLines() {
            Graphics2D graphics2D = (Graphics2D) getGraphics();
            for (Line2D line: lines){
                graphics2D.draw(line);
            }
        }

        private void drawC(capacitor capacitor1) {
            int x1, y1, x2, y2;
            int node1Num = Integer.parseInt(capacitor1.n1.name) - 1;
            int node2Num = Integer.parseInt(capacitor1.n2.name) - 1;
            if (node1Num == -1)
                connectToEarth(node2Num, "CapacitorV.png", capacitor1.NameC);
            else if (node2Num == -1)
                connectToEarth(node1Num, "CapacitorV.png", capacitor1.NameC);
            else {
                x1 = node1Num % 6;
                y1 = 4 - node1Num / 6;
                x2 = node2Num % 6;
                y2 = 4 - node2Num / 6;
                if (x1 == x2) {
                    drawVertical(x1, y1, x2, y2, "CapacitorV.png", capacitor1.NameC);
                } else if (y1 == y2) {
                    drawHorizontal(x1, y1, x2, y2, "CapacitorH.png", capacitor1.NameC);
                }
            }
        }
        private void drawVS(voltageSource voltageSource1) {
            int x1, y1, x2, y2;
            int node1Num = Integer.parseInt(voltageSource1.n1.name) - 1;
            int node2Num = Integer.parseInt(voltageSource1.n2.name) - 1;
            if (node1Num == -1)
                connectToEarth(node2Num, "VoltageU.png", voltageSource1.NameV);
            else if (node2Num == -1)
                connectToEarth(node1Num, "VoltageD.png", voltageSource1.NameV);
            else {
                x1 = node1Num % 6;
                y1 = 4 - node1Num / 6;
                x2 = node2Num % 6;
                y2 = 4 - node2Num / 6;
                if (x1 == x2) {
                    if (y1 > y2)
                        drawVertical(x1, y1, x2, y2, "VoltageU.png", voltageSource1.NameV);
                    else
                        drawVertical(x1, y1, x2, y2, "VoltageD.png", voltageSource1.NameV);
                } else if (y1 == y2) {
                    if (x1 > x2)
                        drawHorizontal(x1, y1, x2, y2, "VoltageL.png", voltageSource1.NameV);
                    else
                        drawHorizontal(x1, y1, x2, y2, "VoltageR.png", voltageSource1.NameV);
                }
            }
        }

        private void drawL(inductor inductor1) {
            int x1, y1, x2, y2;
            int node1Num = Integer.parseInt(inductor1.n1.name) - 1;
            int node2Num = Integer.parseInt(inductor1.n2.name) - 1;
            if (node1Num == -1)
                connectToEarth(node2Num, "InductorV.png", inductor1.NameL);
            else if (node2Num == -1)
                connectToEarth(node1Num, "InductorV.png", inductor1.NameL);
            else {
                x1 = node1Num % 6;
                y1 = 4 - node1Num / 6;
                x2 = node2Num % 6;
                y2 = 4 - node2Num / 6;
                if (x1 == x2) {
                    drawVertical(x1, y1, x2, y2, "InductorV.png", inductor1.NameL);
                } else if (y1 == y2) {
                    drawHorizontal(x1, y1, x2, y2, "InductorH.png", inductor1.NameL);
                }
            }
        }

        private void drawCS(currentSource currentSource1) {
            int x1, y1, x2, y2;
            int node1Num = Integer.parseInt(currentSource1.n1.name) - 1;
            int node2Num = Integer.parseInt(currentSource1.n2.name) - 1;
            if (node1Num == -1)
                connectToEarth(node2Num, "CurrentU.png", currentSource1.NameI);
            else if (node2Num == -1)
                connectToEarth(node1Num, "CurrentD.png", currentSource1.NameI);
            else {
                x1 = node1Num % 6;
                y1 = 4 - node1Num / 6;
                x2 = node2Num % 6;
                y2 = 4 - node2Num / 6;
                if (x1 == x2) {
                    if (y1 > y2)
                        drawVertical(x1, y1, x2, y2, "CurrentU.png", currentSource1.NameI);
                    else
                        drawVertical(x1, y1, x2, y2, "CurrentD.png", currentSource1.NameI);
                } else if (y1 == y2) {
                    if (x1 > x2)
                        drawHorizontal(x1, y1, x2, y2, "CurrentL.png", currentSource1.NameI);
                    else
                        drawHorizontal(x1, y1, x2, y2, "CurrentR.png", currentSource1.NameI);
                }
            }
        }
        private void drawResistor(resistor resistor1) {
            int x1, y1, x2, y2;
            int node1Num = Integer.parseInt(resistor1.n1.name) - 1;
            int node2Num = Integer.parseInt(resistor1.n2.name) - 1;
            if (node1Num == -1)
                connectToEarth(node2Num, "ResistorV.png", resistor1.NameR);
            else if (node2Num == -1)
                connectToEarth(node1Num, "ResistorV.png", resistor1.NameR);
            else {
                x1 = node1Num % 6;
                y1 = 4 - node1Num / 6;
                x2 = node2Num % 6;
                y2 = 4 - node2Num / 6;
                System.out.println(x1 + " " + y1 + " " + x2 + " " + y2);
                if (x1 == x2) {
                    drawVertical(x1, y1, x2, y2, "ResistorV.png", resistor1.NameR);
                } else if (y1 == y2) {
                    drawHorizontal(x1, y1, x2, y2, "ResistorH.png", resistor1.NameR);
                }
            }
        }
        private void connectToEarth(int nodeNum, String address, String name) {
            int x1, x2, y1, y2;
            x1 = nodeNum % 6;
            y1 = 5;
            y2 = 4 - nodeNum / 6;
            drawVertical(x1, y1, x1, y2, address, name);
        }

        private void drawHorizontal(int x1, int y1, int x2, int y2, String address, String name) {
            if (x1 > x2){
                int temp = x2;
                x2 = x1;
                x1 = temp;
            }
            String string = x1 + "" + y1 + "" + x2 + "" + y2;
            int count = 0;
            if (paralel.keySet().contains(string))
                count = paralel.get(string);
            BufferedImage img = null;
            try {
                img = ImageIO.read(new File(address));
            } catch (IOException e) {
                e.printStackTrace();
            }
            JLabel jLabel = new JLabel();
            if (address.contains("Inductor"))
                jLabel.setBounds((int)((double)(x1 + x2) / 2.0 * (1200 / 6)) + 35, y1 * 800 / 6 + count * 40 + 12, 70, 25);
            else if (address.contains("Current"))
                jLabel.setBounds((int)((double)(x1 + x2) / 2.0 * (1200 / 6)) + 35, y1 * 800 / 6 + count * 40 + 10, 70, 32);
            else
                jLabel.setBounds((int)((double)(x1 + x2) / 2.0 * (1200 / 6)) + 35, y1 * 800 / 6 + count * 40, 70, 50);
            Image image;
            image = img.getScaledInstance(jLabel.getWidth(), jLabel.getHeight(), Image.SCALE_SMOOTH);
            ImageIcon imageIcon = new ImageIcon(image);
            jLabel.setIcon(imageIcon);
            add(jLabel);
            super.paint(getGraphics());

            int additionalX = 0, additionalY = 0;
            if (address.contains("Capacitor"))
            {
                additionalX += 3;
                additionalY -= 2;
            }
            if (address.contains("Inductor")){
                additionalY -= 4;
            }
            if (address.contains("Current"))
            {
                additionalX += 3;
                additionalY += 9;
            }
            if (address.contains("Voltage")){
                additionalX -= 10;
            }
            JLabel jLabel1 = new JLabel(name);
            jLabel1.setBounds((int)((double)(x1 + x2) / 2.0 * (1200 / 6)) + 60 + additionalX, y1 * 800 / 6 + count * 40 - 15 + additionalY, 70, 50);

            add(jLabel1);
            super.paint(getGraphics());

            Line2D line1 = new Line2D.Float(x1 * 1200 / 6 + 100, y1 * 800 / 6 + 80, (int)((double)(x1 + x2) / 2.0 * (1200 / 6)) + 45, y1 * 800 / 6 + count * 40 + 55+8);
            Line2D line2 = new Line2D.Float((int)((double)(x1 + x2) / 2.0 * (1200 / 6)) + 110, y1 * 800 / 6 + count * 40 + 55+8, x2 * 1200 / 6 + 100, y1 * 800 / 6 + 80);
            lines.add(line1);
            lines.add(line2);
            paralel.put(string, count + 1);

        }

        private void drawVertical(int x1, int y1, int x2, int y2, String address, String name) {
            if (y1 > y2) {
                int temp = y1;
                y1 = y2;
                y2 = temp;
            }
            String string = x1 + "" + y1 + "" + x2 + "" + y2;
            int count = 0;
            if (paralel.keySet().contains(string))
                count = paralel.get(string);
            BufferedImage img = null;
            try {
                img = ImageIO.read(new File(address));
            } catch (IOException e) {
                e.printStackTrace();
            }

            JLabel jLabel = new JLabel();
            if (address.contains("Inductor"))
                jLabel.setBounds(x1 * 1200 / 6 + count * 40 + 22, (int)((double)(y1 + y2) / 2.0 * (800 / 6)) - 5, 25, 45);
            else if (address.contains("Current"))
                jLabel.setBounds(x1 * 1200 / 6 + count * 40 + 15, (int)((double)(y1 + y2) / 2.0 * (800 / 6)) - 5, 40, 50);
            else if (address.contains("Voltage"))
                jLabel.setBounds(x1 * 1200 / 6 + count * 40 + 17, (int)((double)(y1 + y2) / 2.0 * (800 / 6)) - 5, 35, 50);
            else
                jLabel.setBounds(x1 * 1200 / 6 + count * 40, (int)((double)(y1 + y2) / 2.0 * (800 / 6)) - 5, 70, 50);
            Image image = img.getScaledInstance(jLabel.getWidth(), jLabel.getHeight(), Image.SCALE_SMOOTH);
            ImageIcon imageIcon = new ImageIcon(image);
            jLabel.setIcon(imageIcon);
            add(jLabel);

            super.paint(getGraphics());

            int additionalX = 0, additionalY = 0;
            if (address.contains("Capacitor"))
            {
                additionalX += 15;
            }
            if (address.contains("Current"))
            {
                additionalX += 13;
            }
            if (address.contains("Voltage")){
                additionalX += 7;
                additionalY -= 11;
            }
            JLabel jLabel1 = new JLabel(name);
            jLabel1.setBounds(x1 * 1200 / 6 + count * 40 + 8 + additionalX, (int)((double)(y1 + y2) / 2.0 * (800 / 6)) - 5 + additionalY, 70, 50);

            add(jLabel1);
            super.paint(getGraphics());
            Line2D line1 = new Line2D.Float(x1 * 1200 / 6 + 100, y1 * 800 / 6 + 80, x1 * 1200 / 6 + count * 40 + 42, (int)((double)(y1 + y2) / 2.0 * (800 / 6)) + 27+8);
            lines.add(line1);
            Line2D line2 = new Line2D.Float(x1 * 1200 / 6 + count * 40 + 42, (int)((double)(y1 + y2) / 2.0 * (800 / 6)) + 71+8, x2 * 1200 / 6 + 100, y2 * 800 / 6 + 80);
            lines.add(line2);

            paralel.put(string, count + 1);
        }

    }

    public static void main(String[] args) {

        node n;
        resistor r;
        currentSource cs;
        voltageSource vs;
        capacitor c;
        inductor l;
        CCCS cccs;
        VCCS vccs;
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
                l.outputCurrent.add(0.0);
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
                vs.outputCurrent.add(0.0);
                VS.add(vs);
                input.add("V");
            }
            else if (s.charAt(0) == 'C') {
                c = new capacitor(s);
                c.outputCurrent.add(0.0);
                C.add(c);
                input.add("C");
            }
            else if (s.charAt(0) == 'G') {
                vccs = new VCCS(s);
                G.add(vccs);
                input.add("G");
            }
            else if (s.charAt(0) == 'F') {
                cccs = new CCCS(s);
                F.add(cccs);
                input.add("F");
            }
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
                calcNodeVolts(dT, dV, dI, i);
                calcBranchCurrents(dT, dV, dI, i);
            }
        }

        chapOutput();

        graphProject GraphProject = new graphProject();

    }
}
