
import javafx.stage.FileChooser;

import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.filechooser.FileSystemView;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.Line2D;
import java.awt.image.BufferedImage;
import java.awt.image.ImageObserver;
import java.awt.image.ImageProducer;
import java.io.File;
import java.io.FileWriter;
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
    static ArrayList<VCVS> E=new ArrayList<VCVS>(0);
    static ArrayList<CCVS> H=new ArrayList<CCVS>(0);
    static ArrayList<diode> D=new ArrayList<diode>(0);


    static File file = new File("test1.txt");


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
        else if (s.charAt(0) == 'E') {
            for (int i=0;i<E.size();i++){
                if(E.get(i).NameE.equals(s)){
                    return E.get(i).outputCurrent;
                }
            }
        }
        else if (s.charAt(0) == 'H') {
            for (int i=0;i<H.size();i++){
                if(H.get(i).NameH.equals(s)){
                    return H.get(i).outputCurrent;
                }
            }
        }
        else if (s.charAt(0) == 'D') {
            for (int i=0;i<D.size();i++){
                if(D.get(i).NameD.equals(s)){
                    return D.get(i).outputCurrent;
                }
            }
        }
        s+=0/0;
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
        public int searchNodeInUnion(String s){
            for(int i=0;i<n.size();i++){
                if(n.get(i).name.equals(s)){
                    return i;
                }
            }
            return -1;
        }
        public int addInputUnionVolts(int i, int n_Place, int iteration, double dT){
            if(i<n.size()){
                for(int j=0, k=0;j<VS.size();j++){
                    if(VS.get(j).n1.name.equals(n.get(n_Place).name) && iteration+2>VS.get(j).n2.outputVolt.size()){
                        VS.get(j).n2.volt=VS.get(j).V+VS.get(j).A*Math.sin(VS.get(j).p+VS.get(j).w*(iteration)*dT)+VS.get(j).n1.volt;
                        VS.get(j).n2.outputVolt.add(VS.get(j).n2.volt);
                        k=searchNodeInUnion(VS.get(j).n2.name);
                        i=addInputUnionVolts(i+1, k, iteration, dT);
                    }
                    else if(VS.get(j).n2.name.equals(n.get(n_Place).name) && iteration+2>VS.get(j).n1.outputVolt.size()){
                        VS.get(j).n1.volt=-VS.get(j).V-VS.get(j).A*Math.sin(VS.get(j).p+VS.get(j).w*(iteration)*dT)+VS.get(j).n2.volt;
                        VS.get(j).n1.outputVolt.add(VS.get(j).n1.volt);
                        k=searchNodeInUnion(VS.get(j).n1.name);
                        i=addInputUnionVolts(i+1, k, iteration, dT);
                    }


                }
                for(int j=0, k=0;j<E.size();j++){
                    if(E.get(j).n1.name.equals(n.get(n_Place).name) && iteration+2>E.get(j).n2.outputVolt.size()){
                        E.get(j).updateV(iteration);
                        E.get(j).n2.volt=E.get(j).V+E.get(j).n1.volt;
                        E.get(j).n2.outputVolt.add(E.get(j).n2.volt);
                        k=searchNodeInUnion(E.get(j).n2.name);
                        i=addInputUnionVolts(i+1, k,iteration, dT);
                    }
                    else if(E.get(j).n2.name.equals(n.get(n_Place).name) && iteration+2>E.get(j).n1.outputVolt.size()){
                        E.get(j).updateV(iteration);
                        E.get(j).n1.volt=E.get(j).V+E.get(j).n2.volt;
                        E.get(j).n1.outputVolt.add(E.get(j).n1.volt);
                        k=searchNodeInUnion(E.get(j).n1.name);
                        i=addInputUnionVolts(i+1, k,iteration, dT);
                    }
                }
                for(int j=0, k=0;j<H.size();j++){
                    if(H.get(j).n1.name.equals(n.get(n_Place).name) && iteration+2>H.get(j).n2.outputVolt.size()){
                        H.get(j).updateV(iteration);
                        H.get(j).n2.volt=H.get(j).V+H.get(j).n1.volt;
                        H.get(j).n2.outputVolt.add(H.get(j).n2.volt);
                        k=searchNodeInUnion(H.get(j).n2.name);
                        i=addInputUnionVolts(i+1, k, iteration, dT);
                    }
                    else if(H.get(j).n2.name.equals(n.get(n_Place).name) && iteration+2>H.get(j).n1.outputVolt.size()){
                        H.get(j).updateV(iteration);
                        H.get(j).n1.volt=H.get(j).V+H.get(j).n2.volt;
                        H.get(j).n1.outputVolt.add(H.get(j).n1.volt);
                        k=searchNodeInUnion(H.get(j).n1.name);
                        i=addInputUnionVolts(i+1, k, iteration, dT);
                    }
                }
            }
            return i;
        }
        public int addInputUnionVolts0(int i, int n_Place, int iteration, double dT){
            if(i<n.size()){
                for(int j=0, k=0;j<VS.size();j++){
                    if(VS.get(j).n1.name.equals(n.get(n_Place).name) && iteration+1>VS.get(j).n2.outputVolt.size()){
                        VS.get(j).n2.volt=VS.get(j).V+VS.get(j).A*Math.sin(VS.get(j).p+VS.get(j).w*(iteration)*dT)+VS.get(j).n1.volt;
                        VS.get(j).n2.outputVolt.add(VS.get(j).n2.volt);
                        k=searchNodeInUnion(VS.get(j).n2.name);
                        i=addInputUnionVolts0(i+1, k, iteration, dT);
                    }
                    else if(VS.get(j).n2.name.equals(n.get(n_Place).name) && iteration+1>VS.get(j).n1.outputVolt.size()){
                        VS.get(j).n1.volt=-VS.get(j).V-VS.get(j).A*Math.sin(VS.get(j).p+VS.get(j).w*(iteration)*dT)+VS.get(j).n2.volt;
                        VS.get(j).n1.outputVolt.add(VS.get(j).n1.volt);
                        k=searchNodeInUnion(VS.get(j).n1.name);
                        i=addInputUnionVolts0(i+1, k, iteration, dT);
                    }


                }
                for(int j=0, k=0;j<E.size();j++){
                    if(E.get(j).n1.name.equals(n.get(n_Place).name) && iteration+1>E.get(j).n2.outputVolt.size()){
                        E.get(j).updateV(iteration);
                        E.get(j).n2.volt=E.get(j).V+E.get(j).n1.volt;
                        E.get(j).n2.outputVolt.add(E.get(j).n2.volt);
                        k=searchNodeInUnion(E.get(j).n2.name);
                        i=addInputUnionVolts0(i+1, k,iteration, dT);
                    }
                    else if(E.get(j).n2.name.equals(n.get(n_Place).name) && iteration+1>E.get(j).n1.outputVolt.size()){
                        E.get(j).updateV(iteration);
                        E.get(j).n1.volt=E.get(j).V+E.get(j).n2.volt;
                        E.get(j).n1.outputVolt.add(E.get(j).n1.volt);
                        k=searchNodeInUnion(E.get(j).n1.name);
                        i=addInputUnionVolts0(i+1, k,iteration, dT);
                    }
                }
                for(int j=0, k=0;j<H.size();j++){
                    if(H.get(j).n1.name.equals(n.get(n_Place).name) && iteration+1>H.get(j).n2.outputVolt.size()){
                        H.get(j).updateV(iteration);
                        H.get(j).n2.volt=H.get(j).V+H.get(j).n1.volt;
                        H.get(j).n2.outputVolt.add(H.get(j).n2.volt);
                        k=searchNodeInUnion(H.get(j).n2.name);
                        i=addInputUnionVolts0(i+1, k, iteration, dT);
                    }
                    else if(H.get(j).n2.name.equals(n.get(n_Place).name) && iteration+1>H.get(j).n1.outputVolt.size()){
                        H.get(j).updateV(iteration);
                        H.get(j).n1.volt=H.get(j).V+H.get(j).n2.volt;
                        H.get(j).n1.outputVolt.add(H.get(j).n1.volt);
                        k=searchNodeInUnion(H.get(j).n1.name);
                        i=addInputUnionVolts0(i+1, k, iteration, dT);
                    }
                }
            }
            return i;
        }
    }
    abstract public static class branch{
        ArrayList<Double> outputCurrent=new ArrayList<Double>(0);
        node n1, n2;
        double I=0;
        public void output(FileWriter fileWriter){
            try {
                for (int i = 0; i < n1.outputVolt.size(); i++) {
                    fileWriter.write(" " + (n1.outputVolt.get(i) - n2.outputVolt.get(i)) + "|" + outputCurrent.get(i) + "|" + (outputCurrent.get(i) * (n1.outputVolt.get(i) - n2.outputVolt.get(i))));
                }
            }
            catch (Exception e){

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
            if(n.equal(n1)!=1) {
                if (i == -1) {
                    n.union = n1.union;
                    U.get(j).n.add(n);
                    N.add(n);
                    n2 = n;
                }
                else {
                    int k = searchUnion(N.get(i).union);
                    for (int l = 0; l < U.get(k).n.size(); l++) {
                        U.get(k).n.get(l).union = U.get(j).union;
                        U.get(j).n.add(U.get(k).n.get(l));
                    }
                    U.remove(k);
                    n2 = N.get(i);
                }
            }
            else {
                n2=n1;
            }
            s=s.substring(s.indexOf(" ")+1);
            s=s.replaceAll("k", "000");
            s=s.replaceAll("M", "000000");
            s=s.replaceAll("G", "000000000");
            if(s.substring(0, s.indexOf(" ")).indexOf("m")!=-1){
                V=0.001;
                s=s.replaceAll("m", "");
            }
            else if(s.substring(0, s.indexOf(" ")).indexOf("u")!=-1){
                V=0.000001;
                s=s.replaceAll("u", "");
            }
            else if(s.substring(0, s.indexOf(" ")).indexOf("n")!=-1){
                V=0.000000001;
                s=s.replaceAll("n", "");
            }
            else if(s.substring(0, s.indexOf(" ")).indexOf("p")!=-1){
                V=0.000000000001;
                s=s.replaceAll("p", "");
            }
            else {
                V=1;
            }
            V*=Double.parseDouble(s.substring(0, s.indexOf(" ")));
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
            A=Double.parseDouble(s.substring(0,s.indexOf(" ")));
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
            w=2*Math.PI*Double.parseDouble(s.substring(0,s.indexOf(" ")));
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
            p=Double.parseDouble(s);
            V+=A*Math.sin(p);
        }
        public double calcOutputCurrent(double dT, double dV, double dI, int iteration){
            double I=0;
            for(int k=0; k<R.size();k++){
                if(R.get(k).n1.name.equals(n1.name)){
                    I+=(R.get(k).n2.outputVolt.get(iteration)-R.get(k).n1.outputVolt.get(iteration))/R.get(k).R;
                }
                else if(R.get(k).n2.name.equals(n1.name)){
                    I-=(R.get(k).n2.outputVolt.get(iteration)-R.get(k).n1.outputVolt.get(iteration))/R.get(k).R;
                }
            }
            for(int k=0; k<CS.size();k++){
                if(CS.get(k).n1.name.equals(n1.name)){
                    I += CS.get(k).outputCurrent.get(iteration);
                }
                else if(CS.get(k).n2.name.equals(n1.name)){
                    I-=CS.get(k).outputCurrent.get(iteration);
                }
            }
            for(int k=0; k<C.size();k++){
                if(C.get(k).n1.name.equals(n1.name)) {
                    I-= C.get(k).C * (C.get(k).n1.moshtaghVolt(dT,iteration) - C.get(k).n2.moshtaghVolt(dT,iteration));
                }
                else if(C.get(k).n2.name.equals(n1.name)){
                    I += C.get(k).C * (C.get(k).n1.moshtaghVolt(dT, iteration) - C.get(k).n2.moshtaghVolt(dT, iteration));
                }
            }
            for(int k=0; k<L.size();k++){
                if(L.get(k).n1.name.equals(n1.name)){
                    I-=L.get(k).outputCurrent.get(iteration);
                }
                else if(L.get(k).n2.name.equals(n1.name)){
                    I+=L.get(k).outputCurrent.get(iteration);
                }
            }
            for(int k=0; k<G.size();k++){
                if(G.get(k).n1.name.equals(n1.name)){
                    I+=G.get(k).outputCurrent.get(iteration);
                }
                else if(G.get(k).n2.name.equals(n1.name)){
                    I-=G.get(k).outputCurrent.get(iteration);
                }
            }
            for(int k=0; k<F.size();k++){
                if(F.get(k).n1.name.equals(n1.name)){
                    I+=F.get(k).outputCurrent.get(iteration);
                }
                else if(F.get(k).n2.name.equals(n1.name)){
                    I-=F.get(k).outputCurrent.get(iteration);
                }
            }
            for(int k=0; k<VS.size();k++){
                if(!VS.get(k).NameV.equals(NameV)) {
                    if (VS.get(k).n1.name.equals(n1.name)) {
                        I -= VS.get(k).outputCurrent.get(iteration);
                    }
                    else if (VS.get(k).n2.name.equals(n1.name)) {
                        I += VS.get(k).outputCurrent.get(iteration);
                    }
                }
            }
            for(int k=0; k<E.size();k++) {
                if (E.get(k).n1.name.equals(n1.name)) {
                    I -= E.get(k).outputCurrent.get(iteration);
                } else if (E.get(k).n2.name.equals(n1.name)) {
                    I += E.get(k).outputCurrent.get(iteration);
                }
            }
            for(int k=0; k<H.size();k++) {
                if (H.get(k).n1.name.equals(n1.name)) {
                    I -= H.get(k).outputCurrent.get(iteration);
                } else if (H.get(k).n2.name.equals(n1.name)) {
                    I += H.get(k).outputCurrent.get(iteration);
                }
            }
            for(int k=0; k<D.size();k++){
                if(D.get(k).n1.name.equals(n1.name)){
                    I-=D.get(k).outputCurrent.get(iteration);
                }
                else if(D.get(k).n2.name.equals(n1.name)){
                    I+=D.get(k).outputCurrent.get(iteration);
                }
            }
            return I;
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
                s+=0/0;
            }
            else{
                n3=N.get(i);
            }
            s=s.substring(s.indexOf(" ")+1);
            n=new node(s.substring(0, s.indexOf(" ")));
            i=searchNode(s.substring(0,s.indexOf(" ")));
            if(i==-1){
                s+=0/0;
            }
            else{
                n4=N.get(i);
            }
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
            I=0;
            I*=a;
            outputCurrent.add(I);
        }
    }
    public static class CCCS extends branch {
        String NameF;
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
    public static class VCVS extends branch {
        String NameE;
        node n3, n4;
        double a, V;
        VCVS(String s){
            NameE=s.substring(0,s.indexOf(" "));
            s=s.substring(s.indexOf(" ")+1);
            node n=new node(s.substring(0, s.indexOf(" ")));
            nodes nodz=new nodes();
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
            if(n.equal(n1)!=1) {
                if (i == -1) {
                    n.union = n1.union;
                    U.get(j).n.add(n);
                    N.add(n);
                    n2 = n;
                }
                else {
                    int k = searchUnion(N.get(i).union);
                    for (int l = 0; l < U.get(k).n.size(); l++) {
                        U.get(k).n.get(l).union = U.get(j).union;
                        U.get(j).n.add(U.get(k).n.get(l));
                    }
                    U.remove(k);
                    n2 = N.get(i);
                }
            }
            else {
                n2=n1;
            }
            s=s.substring(s.indexOf(" ")+1);
            n=new node(s.substring(0, s.indexOf(" ")));
            i=searchNode(s.substring(0,s.indexOf(" ")));
            if(i==-1){
                s+=0/0;
            }
            else{
                n3=N.get(i);
            }
            s=s.substring(s.indexOf(" ")+1);
            n=new node(s.substring(0, s.indexOf(" ")));
            i=searchNode(s.substring(0,s.indexOf(" ")));
            if(i==-1){
                s+=0/0;
            }
            else{
                n4=N.get(i);
            }
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
            V=0;
            V*=a;
            outputCurrent.add(0.0);
        }
        public double calcOutputCurrent(double dT, double dV, double dI, int iteration){
            double I=0;
            for(int k=0; k<R.size();k++){
                if(R.get(k).n1.name.equals(n1.name)){
                    I+=(R.get(k).n2.outputVolt.get(iteration)-R.get(k).n1.outputVolt.get(iteration))/R.get(k).R;
                }
                else if(R.get(k).n2.name.equals(n1.name)){
                    I-=(R.get(k).n2.outputVolt.get(iteration)-R.get(k).n1.outputVolt.get(iteration))/R.get(k).R;
                }
            }
            for(int k=0; k<CS.size();k++){
                if(CS.get(k).n1.name.equals(n1.name)){
                    I += CS.get(k).outputCurrent.get(iteration);
                }
                else if(CS.get(k).n2.name.equals(n1.name)){
                    I-=CS.get(k).outputCurrent.get(iteration);
                }
            }
            for(int k=0; k<C.size();k++){
                if(C.get(k).n1.name.equals(n1.name)) {
                    I-= C.get(k).C * (C.get(k).n1.moshtaghVolt(dT,iteration) - C.get(k).n2.moshtaghVolt(dT,iteration));
                }
                else if(C.get(k).n2.name.equals(n1.name)){
                    I += C.get(k).C * (C.get(k).n1.moshtaghVolt(dT, iteration) - C.get(k).n2.moshtaghVolt(dT, iteration));
                }
            }
            for(int k=0; k<L.size();k++){
                if(L.get(k).n1.name.equals(n1.name)){
                    I-=L.get(k).outputCurrent.get(iteration);
                }
                else if(L.get(k).n2.name.equals(n1.name)){
                    I+=L.get(k).outputCurrent.get(iteration);
                }
            }
            for(int k=0; k<G.size();k++){
                if(G.get(k).n1.name.equals(n1.name)){
                    I+=G.get(k).outputCurrent.get(iteration);
                }
                else if(G.get(k).n2.name.equals(n1.name)){
                    I-=G.get(k).outputCurrent.get(iteration);
                }
            }
            for(int k=0; k<F.size();k++){
                if(F.get(k).n1.name.equals(n1.name)){
                    I+=F.get(k).outputCurrent.get(iteration);
                }
                else if(F.get(k).n2.name.equals(n1.name)){
                    I-=F.get(k).outputCurrent.get(iteration);
                }
            }
            for(int k=0; k<VS.size();k++) {
                if (VS.get(k).n1.name.equals(n1.name)) {
                    I -= VS.get(k).outputCurrent.get(iteration);
                }
                else if (VS.get(k).n2.name.equals(n1.name)) {
                    I += VS.get(k).outputCurrent.get(iteration);
                }
            }
            for(int k=0; k<E.size();k++) {
                if(!E.get(k).NameE.equals(NameE)) {
                    if (E.get(k).n1.name.equals(n1.name)) {
                        I -= E.get(k).outputCurrent.get(iteration);
                    } else if (E.get(k).n2.name.equals(n1.name)) {
                        I += E.get(k).outputCurrent.get(iteration);
                    }
                }
            }
            for(int k=0; k<H.size();k++) {
                if (H.get(k).n1.name.equals(n1.name)) {
                    I -= H.get(k).outputCurrent.get(iteration);
                } else if (H.get(k).n2.name.equals(n1.name)) {
                    I += H.get(k).outputCurrent.get(iteration);
                }
            }
            for(int k=0; k<D.size();k++){
                if(D.get(k).n1.name.equals(n1.name)){
                    I-=D.get(k).outputCurrent.get(iteration);
                }
                else if(D.get(k).n2.name.equals(n1.name)){
                    I+=D.get(k).outputCurrent.get(iteration);
                }
            }
            return I;
        }
        public void updateV(int iteration){
            if(n3.outputVolt.size()>0&&n4.outputVolt.size()>0)
                V=n3.outputVolt.get(iteration)-n4.outputVolt.get(iteration);
            else
                V=0;
            V*=a;
        }
    }
    public static class CCVS extends branch {
        String NameH;
        ArrayList<Double> inputCurrent=new ArrayList<Double>(0);
        double a, V=0;
        CCVS(String s){
            NameH=s.substring(0,s.indexOf(" "));
            s=s.substring(s.indexOf(" ")+1);
            node n=new node(s.substring(0, s.indexOf(" ")));
            nodes nodz=new nodes();
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
            if(n.equal(n1)!=1) {
                if (i == -1) {
                    n.union = n1.union;
                    U.get(j).n.add(n);
                    N.add(n);
                    n2 = n;
                }
                else {
                    int k = searchUnion(N.get(i).union);
                    for (int l = 0; l < U.get(k).n.size(); l++) {
                        U.get(k).n.get(l).union = U.get(j).union;
                        U.get(j).n.add(U.get(k).n.get(l));
                    }
                    U.remove(k);
                    n2 = N.get(i);
                }
            }
            else {
                n2=n1;
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
            V=inputCurrent.get(0);
            V*=a;
            outputCurrent.add(0.0);
        }
        public double calcOutputCurrent(double dT, double dV, double dI, int iteration){
            double I=0;
            for(int k=0; k<R.size();k++){
                if(R.get(k).n1.name.equals(n1.name)){
                    I+=(R.get(k).n2.outputVolt.get(iteration)-R.get(k).n1.outputVolt.get(iteration))/R.get(k).R;
                }
                else if(R.get(k).n2.name.equals(n1.name)){
                    I-=(R.get(k).n2.outputVolt.get(iteration)-R.get(k).n1.outputVolt.get(iteration))/R.get(k).R;
                }
            }
            for(int k=0; k<CS.size();k++){
                if(CS.get(k).n1.name.equals(n1.name)){
                    I += CS.get(k).outputCurrent.get(iteration);
                }
                else if(CS.get(k).n2.name.equals(n1.name)){
                    I-=CS.get(k).outputCurrent.get(iteration);
                }
            }
            for(int k=0; k<C.size();k++){
                if(C.get(k).n1.name.equals(n1.name)) {
                    I-= C.get(k).C * (C.get(k).n1.moshtaghVolt(dT,iteration) - C.get(k).n2.moshtaghVolt(dT,iteration));
                }
                else if(C.get(k).n2.name.equals(n1.name)){
                    I += C.get(k).C * (C.get(k).n1.moshtaghVolt(dT, iteration) - C.get(k).n2.moshtaghVolt(dT, iteration));
                }
            }
            for(int k=0; k<L.size();k++){
                if(L.get(k).n1.name.equals(n1.name)){
                    I-=L.get(k).outputCurrent.get(iteration);
                }
                else if(L.get(k).n2.name.equals(n1.name)){
                    I+=L.get(k).outputCurrent.get(iteration);
                }
            }
            for(int k=0; k<G.size();k++){
                if(G.get(k).n1.name.equals(n1.name)){
                    I+=G.get(k).outputCurrent.get(iteration);
                }
                else if(G.get(k).n2.name.equals(n1.name)){
                    I-=G.get(k).outputCurrent.get(iteration);
                }
            }
            for(int k=0; k<F.size();k++){
                if(F.get(k).n1.name.equals(n1.name)){
                    I+=F.get(k).outputCurrent.get(iteration);
                }
                else if(F.get(k).n2.name.equals(n1.name)){
                    I-=F.get(k).outputCurrent.get(iteration);
                }
            }
            for(int k=0; k<VS.size();k++) {
                if (VS.get(k).n1.name.equals(n1.name)) {
                    I -= VS.get(k).outputCurrent.get(iteration);
                } else if (VS.get(k).n2.name.equals(n1.name)) {
                    I += VS.get(k).outputCurrent.get(iteration);
                }
            }
            for(int k=0; k<E.size();k++) {
                if (E.get(k).n1.name.equals(n1.name)) {
                    I -= E.get(k).outputCurrent.get(iteration);
                } else if (E.get(k).n2.name.equals(n1.name)) {
                    I += E.get(k).outputCurrent.get(iteration);
                }
            }
            for(int k=0; k<H.size();k++) {
                if(!H.get(k).NameH.equals(NameH)) {
                    if (H.get(k).n1.name.equals(n1.name)) {
                        I -= H.get(k).outputCurrent.get(iteration);
                    }
                    else if (H.get(k).n2.name.equals(n1.name)) {
                        I += H.get(k).outputCurrent.get(iteration);
                    }
                }
            }
            for(int k=0; k<D.size();k++){
                if(D.get(k).n1.name.equals(n1.name)){
                    I-=D.get(k).outputCurrent.get(iteration);
                }
                else if(D.get(k).n2.name.equals(n1.name)){
                    I+=D.get(k).outputCurrent.get(iteration);
                }
            }
            return I;
        }
        public void updateV(int iteration){
            V=inputCurrent.get(iteration);
            V*=a;
        }
    }
    public static class diode extends branch {
        String NameD;
        diode(String s){
            NameD=s.substring(0,s.indexOf(" "));
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
            n=new node(s);
            i=searchNode(s);
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
        }
        public double calcI(int iteration, double dV){
            double i=0, is=0.00001, vt=0.0026;
            i=is*(Math.exp((n1.outputVolt.get(iteration)-n2.outputVolt.get(iteration)+dV)/vt)-1);
            return i;
        }
    }




    public static void calcNodeVolts(double dT, double dV, double dI, int iteration){
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
                        i1+=CS.get(k).outputCurrent.get(iteration);
                    }
                    else if(CS.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i1-=CS.get(k).outputCurrent.get(iteration);
                    }
                }
                for(int k=0; k<C.size();k++){
                    if(C.get(k).n1.name.equals(U.get(i).n.get(j).name)) {
                        i1 -= C.get(k).C * (C.get(k).n1.moshtaghVolt(dT,iteration) - C.get(k).n2.moshtaghVolt(dT,iteration));
                    }
                    else if(C.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i1 += C.get(k).C * (C.get(k).n1.moshtaghVolt(dT, iteration) - C.get(k).n2.moshtaghVolt(dT, iteration));
                    }
                }
                for(int k=0; k<L.size();k++){
                    if(L.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i1-=L.get(k).outputCurrent.get(iteration);
                    }
                    else if(L.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i1+=L.get(k).outputCurrent.get(iteration);
                    }
                }
                for(int k=0; k<G.size();k++){
                    if(G.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i1+=G.get(k).outputCurrent.get(iteration);
                    }
                    else if(G.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i1-=G.get(k).outputCurrent.get(iteration);
                    }
                }
                for(int k=0; k<F.size();k++){
                    if(F.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i1+=F.get(k).outputCurrent.get(iteration);
                    }
                    else if(F.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i1-=F.get(k).outputCurrent.get(iteration);
                    }
                }
                for(int k=0; k<VS.size();k++){
                    if(VS.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i1+=VS.get(k).outputCurrent.get(iteration);
                    }
                    else if(VS.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i1-=VS.get(k).outputCurrent.get(iteration);
                    }
                }
                for(int k=0; k<E.size();k++){
                    if(E.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i1+=E.get(k).outputCurrent.get(iteration);
                    }
                    else if(E.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i1-=E.get(k).outputCurrent.get(iteration);
                    }
                }
                for(int k=0; k<H.size();k++){
                    if(H.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i1+=H.get(k).outputCurrent.get(iteration);
                    }
                    else if(H.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i1-=H.get(k).outputCurrent.get(iteration);
                    }
                }
                for(int k=0; k<D.size();k++){
                    if(D.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i1-=D.get(k).calcI(iteration,0);
                    }
                    else if(D.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i1+=D.get(k).calcI(iteration, 0);
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
                        i2+=CS.get(k).outputCurrent.get(iteration);
                    }
                    else if(CS.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i2-=CS.get(k).outputCurrent.get(iteration);
                    }
                }
                for(int k=0; k<C.size();k++){
                    if(C.get(k).n1.name.equals(U.get(i).n.get(j).name)) {
                        i2 -= C.get(k).C * (C.get(k).n1.moshtaghVolt(dT, iteration)+dV/dT - C.get(k).n2.moshtaghVolt(dT, iteration));
                    }
                    else if(C.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i2 += C.get(k).C * (C.get(k).n1.moshtaghVolt(dT, iteration) - C.get(k).n2.moshtaghVolt(dT, iteration)-dV/dT);
                    }
                }
                for(int k=0; k<L.size();k++){
                    if(L.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i2-=L.get(k).outputCurrent.get(iteration)+dV*dT/L.get(k).L;
                    }
                    else if(L.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i2+=L.get(k).outputCurrent.get(iteration)-dV*dT/L.get(k).L;
                    }
                }
                for(int k=0; k<G.size();k++){
                    if(G.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i2+=G.get(k).outputCurrent.get(iteration);
                        if(G.get(k).n3.name.equals(U.get(i).n.get(j).name)){
                            i2+=G.get(k).a*dV;
                        }
                        else if(G.get(k).n4.name.equals(U.get(i).n.get(j).name)){
                            i2-=G.get(k).a*dV;
                        }
                    }
                    else if(G.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i2-=G.get(k).outputCurrent.get(iteration);
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
                        i2+=F.get(k).outputCurrent.get(iteration);
                    }
                    else if(F.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i2-=F.get(k).outputCurrent.get(iteration);
                    }
                }
                for(int k=0; k<VS.size();k++){
                    if(VS.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i2+=VS.get(k).outputCurrent.get(iteration);
                    }
                    else if(VS.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i2-=VS.get(k).outputCurrent.get(iteration);
                    }
                }
                for(int k=0; k<E.size();k++){
                    if(E.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i2+=E.get(k).outputCurrent.get(iteration);
                    }
                    else if(E.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i2-=E.get(k).outputCurrent.get(iteration);
                    }
                }
                for(int k=0; k<H.size();k++){
                    if(H.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i2+=H.get(k).outputCurrent.get(iteration);
                    }
                    else if(H.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i2-=H.get(k).outputCurrent.get(iteration);
                    }
                }
                for(int k=0; k<D.size();k++){
                    if(D.get(k).n1.name.equals(U.get(i).n.get(j).name)){
                        i2-=D.get(k).calcI(iteration, dV);
                    }
                    else if(D.get(k).n2.name.equals(U.get(i).n.get(j).name)){
                        i2+=D.get(k).calcI(iteration, -dV);
                    }
                }
            }

            U.get(i).n.get(0).outputVolt.add(U.get(i).n.get(0).volt+((Math.abs(i1)-Math.abs(i2))*dV)/dI);
            U.get(i).addInputUnionVolts(0, 0, iteration, dT);
            i1=0;
            i2=0;
        }

        for(int i=0, j;i<U.size();i++){
            j=U.get(i).n.get(0).outputVolt.size()-1;
            U.get(i).n.get(0).volt=U.get(i).n.get(0).outputVolt.get(j);
            U.get(i).addInputUnionVolts(0, 0, iteration, dT);
        }
    }
    public static void calcBranchCurrents(double dT, double dV, double dI, int iteration){
        double I=0;
        for(int i=0;i<R.size();i++){
            R.get(i).outputCurrent.add((R.get(i).n1.volt-R.get(i).n2.volt)/R.get(i).R);
        }
        for (int i=0;i<CS.size();i++){
            CS.get(i).outputCurrent.add(CS.get(i).I+CS.get(i).A*Math.sin(CS.get(i).p+CS.get(i).w*(iteration)*dT));
        }
        for (int i=0; i<C.size();i++){
            C.get(i).outputCurrent.add(C.get(i).C * (C.get(i).n1.moshtaghVolt(dT, iteration) - C.get(i).n2.moshtaghVolt(dT, iteration)));
        }
        for (int i=0; i<L.size();i++){
            I=L.get(i).outputCurrent.get(iteration);
            L.get(i).outputCurrent.add(I+(dT*(L.get(i).n1.outputVolt.get(iteration)-L.get(i).n2.outputVolt.get(iteration)))/L.get(i).L);
        }
        for (int i=0;i<G.size();i++){
            G.get(i).outputCurrent.add((G.get(i).n3.outputVolt.get(iteration)-G.get(i).n4.outputVolt.get(iteration))*G.get(i).a);
        }
        for (int i=0;i<F.size();i++){
            F.get(i).outputCurrent.add(F.get(i).inputCurrent.get(iteration));
        }
        for (int i=0;i<VS.size();i++){
            I=VS.get(i).calcOutputCurrent(dT, dV, dI, iteration);
            VS.get(i).outputCurrent.add(I);
        }
        for (int i=0;i<E.size();i++){
            I=E.get(i).calcOutputCurrent(dT, dV, dI, iteration);
            E.get(i).outputCurrent.add(I);
        }
        for (int i=0;i<H.size();i++){
            H.get(i).outputCurrent.add(I);
        }
        for (int i=0; i<D.size();i++){
            D.get(i).outputCurrent.add(D.get(i).calcI(iteration, 0));
        }
    }


    public static void chapOutput(FileWriter fileWriter) {
        try {
            int groundIndex=searchNode("0");
            for (int i = 0; i < N.size(); i++) {
                fileWriter.write(""+N.get(i).name + " :");
                for (int j = 0; j < N.get(i).outputVolt.size(); j++) {

                    fileWriter.write(" " + (N.get(i).outputVolt.get(j)-N.get(groundIndex).outputVolt.get(j)));

                }
                fileWriter.write("\n");
            }


            for (int i = 0, iR = 0, iI = 0, iC = 0, iL = 0, iV = 0, iG = 0, iF = 0, iE = 0, iH = 0, iD=0; i < input.size(); i++) {
                if (input.get(i).equals("R")) {
                    fileWriter.write(R.get(iR).NameR + " :");
                    R.get(iR).output(fileWriter);
                    iR++;
                }
                else if (input.get(i).equals("I")) {
                    fileWriter.write(CS.get(iI).NameI + " :");
                    CS.get(iI).output(fileWriter);
                    iI++;
                }
                else if (input.get(i).equals("C")) {
                    fileWriter.write(C.get(iC).NameC + " :");
                    C.get(iC).output(fileWriter);
                    iC++;
                }
                else if (input.get(i).equals("L")) {
                    fileWriter.write(L.get(iL).NameL + " :");
                    L.get(iL).output(fileWriter);
                    iL++;
                }
                else if (input.get(i).equals("G")) {
                    fileWriter.write(G.get(iG).NameG + " :");
                    G.get(iG).output(fileWriter);
                    iG++;
                }
                else if (input.get(i).equals("F")) {
                    fileWriter.write(F.get(iF).NameF + " :");
                    F.get(iF).output(fileWriter);
                    iF++;
                }
                else if (input.get(i).equals("V")) {
                    fileWriter.write(VS.get(iV).NameV + " :");
                    VS.get(iV).output(fileWriter);
                    iV++;
                }
                else if (input.get(i).equals("E")) {
                    fileWriter.write(E.get(iE).NameE + " :");
                    E.get(iE).output(fileWriter);
                    iE++;
                }
                else if (input.get(i).equals("H")) {
                    fileWriter.write(H.get(iH).NameH + " :");
                    H.get(iH).output(fileWriter);
                    iH++;
                }
                else if (input.get(i).equals("D")) {
                    fileWriter.write(D.get(iD).NameD + " :");
                    D.get(iD).output(fileWriter);
                    iD++;
                }
                fileWriter.write("\n");
            }
        }
        catch (Exception e) {
        }
    }
    public static int consoleInput(String s, double dT, double T){
        int i=-1, iteration, error=1;
        double t;
        node n1, n2;
        n1=new node(s.substring(0,s.indexOf(" ")));
        i=searchNode(s.substring(0,s.indexOf(" ")));
        if(i==-1){
            error=0;
        }
        else{
            n1=N.get(i);
        }
        s=s.substring(s.indexOf(" ")+1);
        n2=new node(s.substring(0,s.indexOf(" ")));
        i=searchNode(s.substring(0,s.indexOf(" ")));
        if(i==-1){
            error=0;
        }
        else{
            n2=N.get(i);
        }
        s=s.substring(s.indexOf(" ")+1);
        t=Double.parseDouble(s);
        if(t>T||t<0){
            error=0;
        }
        iteration=(int)(t/dT);
        if(error==0){
            return error;
        }
        System.out.println((n1.outputVolt.get(iteration)-n2.outputVolt.get(iteration)-1));
        return error;
    }


    public static int checkGround(double dT, int iteration){
        for(int i=0;i<VS.size();i++){
            if(VS.get(i).n1.name.equals(VS.get(i).n2.name)){
                if (VS.get(i).V+VS.get(i).A*Math.sin(VS.get(i).p+VS.get(i).w*(iteration)*dT)!=0) {
                    return -1;
                }
            }
        }
        for(int i=0;i<E.size();i++){
            if(E.get(i).n1.name.equals(E.get(i).n2.name)){
                if (E.get(i).V!=0) {
                    return -1;
                }
            }
        }
        for(int i=0;i<H.size();i++){
            if(H.get(i).n1.name.equals(H.get(i).n2.name)){
                if (H.get(i).V!=0) {
                    return -1;
                }
            }
        }
        return 0;
    }
    public static int checkVS(double dT, int iteration){
        double v=0;
        for(int i=0;i<U.size();i++){
            for(int j1=0;j1<U.get(i).n.size();j1++){
                for(int j2=j1+1, k=0; j2<U.get(i).n.size();j2++) {
                    for (int m = 0; m < VS.size(); m++) {
                        if (VS.get(m).n1.name.equals(U.get(i).n.get(j1).name)&&VS.get(m).n2.name.equals(U.get(i).n.get(j2).name)) {
                            if (k > 1 && v != VS.get(i).V + VS.get(i).A * Math.sin(VS.get(i).p + VS.get(i).w * (iteration) * dT)) {
                                return -1;
                            }
                            k++;
                            v = VS.get(i).V + VS.get(i).A * Math.sin(VS.get(i).p + VS.get(i).w * (iteration) * dT);
                        }
                        else if (VS.get(m).n2.name.equals(U.get(i).n.get(j1).name)&&VS.get(m).n1.name.equals(U.get(i).n.get(j2).name)) {
                            if (k > 1 && v != -VS.get(i).V + VS.get(i).A * Math.sin(VS.get(i).p + VS.get(i).w * (iteration) * dT)) {
                                return -1;
                            }
                            k++;
                            v = -VS.get(i).V + VS.get(i).A * Math.sin(VS.get(i).p + VS.get(i).w * (iteration) * dT);
                        }
                    }
                    for (int m = 0; m < E.size(); m++) {
                        if (E.get(m).n1.name.equals(U.get(i).n.get(j1).name)&&E.get(m).n2.name.equals(U.get(i).n.get(j2).name)) {
                            if (k > 1 && v != E.get(i).V) {
                                return -1;
                            }
                            k++;
                            v = E.get(i).V;
                        }
                        else if (E.get(m).n2.name.equals(U.get(i).n.get(j1).name)&&E.get(m).n1.name.equals(U.get(i).n.get(j2).name)) {
                            if (k > 1 && v != -E.get(i).V) {
                                return -1;
                            }
                            k++;
                            v = -E.get(i).V;
                        }
                    }
                    for (int m = 0; m < H.size(); m++) {
                        if (H.get(m).n1.name.equals(U.get(i).n.get(j1).name)&&H.get(m).n2.name.equals(U.get(i).n.get(j2).name)) {
                            if (k > 1 && v != H.get(i).V) {
                                return -1;
                            }
                            k++;
                            v = H.get(i).V;
                        }
                        else if (H.get(m).n2.name.equals(U.get(i).n.get(j1).name)&&H.get(m).n1.name.equals(U.get(i).n.get(j2).name)) {
                            if (k > 1 && v != -H.get(i).V) {
                                return -1;
                            }
                            k++;
                            v = -H.get(i).V;
                        }
                    }
                    k = 0;
                }
            }
        }
        return 0;
    }
    public static int checkCS(double dT, int iteration){
        double I=0;
        for(int i=0;i<N.size();i++){
            for(int k=0; k<R.size();k++){
                if(R.get(k).n1.name.equals(N.get(i).name)) {
                    if (!R.get(k).n2.name.equals(N.get(i).name)) {
                        return 0;
                    }
                }
            }
            for(int k=0; k<C.size();k++){
                if(C.get(k).n1.name.equals(N.get(i).name)) {
                    if (!C.get(k).n2.name.equals(N.get(i).name)) {
                        return 0;
                    }
                }
            }
            for(int k=0; k<L.size();k++){
                if(L.get(k).n1.name.equals(N.get(i).name)) {
                    if (!L.get(k).n2.name.equals(N.get(i).name)) {
                        return 0;
                    }
                }
            }
            for(int k=0; k<VS.size();k++){
                if(VS.get(k).n1.name.equals(N.get(i).name)) {
                    if (!VS.get(k).n2.name.equals(N.get(i).name)) {
                        return 0;
                    }
                }
            }
            for(int k=0; k<E.size();k++){
                if(E.get(k).n1.name.equals(N.get(i).name)) {
                    if (!E.get(k).n2.name.equals(N.get(i).name)) {
                        return 0;
                    }
                }
            }
            for(int k=0; k<H.size();k++){
                if(H.get(k).n1.name.equals(N.get(i).name)) {
                    if (!H.get(k).n2.name.equals(N.get(i).name)) {
                        return 0;
                    }
                }
            }
            for(int j=0; j<CS.size();j++){
                if(CS.get(j).n1.name.equals(N.get(i).name)){
                    I+=CS.get(j).outputCurrent.get(iteration);
                }
                else if(CS.get(j).n2.name.equals(N.get(i).name)){
                    I-=CS.get(j).outputCurrent.get(iteration);
                }
            }
            for(int j=0; j<F.size();j++){
                if(F.get(j).n1.name.equals(N.get(i).name)){
                    I+=F.get(j).outputCurrent.get(iteration);
                }
                else if(F.get(j).n2.name.equals(N.get(i).name)){
                    I-=F.get(j).outputCurrent.get(iteration);
                }
            }
            for(int j=0; j<G.size();j++){
                if(G.get(j).n1.name.equals(N.get(i).name)){
                    I+=G.get(j).outputCurrent.get(iteration);
                }
                else if(G.get(j).n2.name.equals(N.get(i).name)){
                    I-=G.get(j).outputCurrent.get(iteration);
                }
            }
            for(int k=0; k<D.size();k++){
                if(D.get(k).n1.name.equals(N.get(i).name)) {
                    if (!D.get(k).n2.name.equals(N.get(i).name)) {
                        return 0;
                    }
                }
            }
            if(I!=0){
                return -1;
            }
            I=0;
        }
        return 0;
    }


    public static class graphProject extends JFrame{

        ArrayList<Line2D> lines = new ArrayList<>();
        HashMap<String, Integer> paralel = new HashMap<>();
        HashMap<String, Integer> paralel2 = new HashMap<>();

        graphProject() {
            //JLabel label=new JLabel(new ImageIcon("7878.png"));
            //label.setBounds(200,350,180,100);
            //add(label);
            JLabel Error=new JLabel("This element is not available.Please try again.");
            Error.setBounds(500,0,400,50);

            JButton draw=new JButton("DRAW");
            draw.setBounds(0,0,200,100);
            add(draw);

            draw.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    String element =JOptionPane.showInputDialog(draw.getParent(),"Enter element name","Choose element",JOptionPane.PLAIN_MESSAGE);
                    int tayid=0;
                    for(int j=0;j<R.size()&&tayid==0;j++)
                        if(R.get(j).NameR.equals(element)){
                            new chartPainterVoltage(R.get(j));
                            new chartPainterCurrent(R.get(j));
                            new chartPainterPower(R.get(j));
                            tayid=1;
                        }
                    for(int j=0;j<C.size()&&tayid==0;j++)
                        if(C.get(j).NameC.equals(element)){
                            new chartPainterVoltage(C.get(j));
                            new chartPainterCurrent(C.get(j));
                            new chartPainterPower(C.get(j));
                            tayid=1;
                        }
                    for(int j=0;j<L.size()&&tayid==0;j++)
                        if(L.get(j).NameL.equals(element)){
                            new chartPainterVoltage(L.get(j));
                            new chartPainterCurrent(L.get(j));
                            new chartPainterPower(L.get(j));
                            tayid=1;
                        }
                    for(int j=0;j<CS.size()&&tayid==0;j++)
                        if(CS.get(j).NameI.equals(element)){
                            new chartPainterVoltage(CS.get(j));
                            new chartPainterCurrent(CS.get(j));
                            new chartPainterPower(CS.get(j));
                            tayid=1;
                        }
                    for(int j=0;j<VS.size()&&tayid==0;j++)
                        if(VS.get(j).NameV.equals(element)){
                            new chartPainterVoltage(VS.get(j));
                            new chartPainterCurrent(VS.get(j));
                            new chartPainterPower(VS.get(j));
                            tayid=1;
                        }
                    for(int j=0;j<G.size()&&tayid==0;j++)
                        if(G.get(j).NameG.equals(element)){
                            new chartPainterVoltage(G.get(j));
                            new chartPainterCurrent(G.get(j));
                            new chartPainterPower(G.get(j));
                            tayid=1;
                        }
                    for(int j=0;j<F.size()&&tayid==0;j++)
                        if(F.get(j).NameF.equals(element)){
                            new chartPainterVoltage(F.get(j));
                            new chartPainterCurrent(F.get(j));
                            new chartPainterPower(F.get(j));
                            tayid=1;
                        }
                    for(int j=0;j<E.size()&&tayid==0;j++)
                        if(E.get(j).NameE.equals(element)){
                            new chartPainterVoltage(E.get(j));
                            new chartPainterCurrent(E.get(j));
                            new chartPainterPower(E.get(j));
                            tayid=1;
                        }
                    for(int j=0;j<H.size()&&tayid==0;j++)
                        if(H.get(j).NameH.equals(element)){
                            new chartPainterVoltage(H.get(j));
                            new chartPainterCurrent(H.get(j));
                            new chartPainterPower(H.get(j));
                            tayid=1;
                        }
                    for(int j=0;j<D.size()&&tayid==0;j++)
                        if(D.get(j).NameD.equals(element)){
                            new chartPainterVoltage(D.get(j));
                            new chartPainterCurrent(D.get(j));
                            new chartPainterPower(D.get(j));
                            tayid=1;
                        }
                    if(tayid==0)
                        JOptionPane.showMessageDialog(draw.getParent(),"This element is not available.Please try again.","Show Error",JOptionPane.INFORMATION_MESSAGE);



                }
            });
            setTitle("Circuit Graph");
            setSize(2000, 1000);
            setLayout(null);
            setVisible(true);
        }

        @Override
        public void paint(Graphics g) {
            super.paint(g);
            g.setColor(Color.BLACK);
            for (resistor resistor1: R) {
                fillHash2(resistor1.n1.name, resistor1.n2.name);
            }
            for (currentSource currentSource1: CS){
                fillHash2(currentSource1.n1.name, currentSource1.n2.name);
            }
            for (VCCS vccs: G){
                fillHash2(vccs.n1.name, vccs.n2.name);
            }
            for (CCCS cccs: F){
                fillHash2(cccs.n1.name, cccs.n2.name);
            }
            for (VCVS vcvs: E){
                fillHash2(vcvs.n1.name, vcvs.n2.name);
            }
            for (CCVS ccvs: H){
                fillHash2(ccvs.n1.name, ccvs.n2.name);
            }
            for (voltageSource voltageSource1: VS){
                fillHash2(voltageSource1.n1.name, voltageSource1.n2.name);
            }
            for (capacitor capacitor1: C){
                fillHash2(capacitor1.n1.name, capacitor1.n2.name);
            }
            for (inductor inductor1: L){
                fillHash2(inductor1.n1.name, inductor1.n2.name);
            }
            for (diode diode1: D){
                fillHash2(diode1.n1.name, diode1.n2.name);
            }
            for (resistor resistor1: R) {
                drawResistor(resistor1);
            }
            for (currentSource currentSource1: CS){
                drawCS(currentSource1);
            }
            for (VCCS vccs: G){
                drawVCCS(vccs);
            }
            for (CCCS cccs: F){
                drawCCCS(cccs);
            }
            for (voltageSource voltageSource1: VS){
                drawVS(voltageSource1);
            }
            for (VCVS vcvs: E){
                drawVCVS(vcvs);
            }
            for (CCVS ccvs: H){
                drawCCVS(ccvs);
            }
            for (capacitor capacitor1: C){
                drawC(capacitor1);
            }
            for (inductor inductor1: L){
                drawL(inductor1);
            }
            for (diode diode1: D){
                drawD(diode1);
            }
            drawLines();
            for (int j=0;j<N.size();j++){
                if (!N.get(j).name.equals("0"))
                    g.drawString(N.get(j).name,100+(Integer.parseInt(N.get(j).name)-1)%6*200,80+(4-Integer.parseInt(N.get(j).name)/6)*800/6);

            }
        }

        private void fillHash2(String name1, String name2) {
            int x1, y1, x2, y2;
            String string = "";
            int node1Num = Integer.parseInt(name1) - 1;
            int node2Num = Integer.parseInt(name2) - 1;
            if (node1Num == -1) {
                x1 = node2Num % 6;
                y1 = 5;
            }
            else {
                x1 = node1Num % 6;
                y1 = 4 - node1Num / 6;
            }
            if (node2Num == -1){
                x2 = node1Num % 6;
                y2 = 5;
            }
            else {
                x2 = node2Num % 6;
                y2 = 4 - node2Num / 6;
            }if (x1 == x2) {
                if (y1 > y2)
                    string = x2 + "" + y2 + "" + x1 + "" + y1;
                else
                    string = x1 + "" + y1 + "" + x2 + "" + y2;
            } else if (y1 == y2) {
                if (x1 > x2)
                    string = x2 + "" + y2 + "" + x1 + "" + y1;
                else
                    string = x1 + "" + y1 + "" + x2 + "" + y2;
            }
            if (!paralel2.keySet().contains(string))
                paralel2.put(string, 1);
            else
                paralel2.put(string, paralel2.get(string) + 1);
        }

        private void drawLines() {
            int min=-2,max=-2;
            for(int j=0;j<6;j++){
                for (String string:paralel2.keySet()){
                    if(Pattern.matches("\\d\\d"+j+"5",string)){
                        max=j;
                        if (min==-2){
                            min=j;
                        }
                    }
                }

            }
            for (int i = min; i < max; i++){
                Line2D line1 = new Line2D.Float(i * 1200 / 6 + 100, 5 * 800 / 6 + 80, (i+1) * 1200 / 6 + 100, 5 * (800 / 6) + 80);
                lines.add(line1);
            }
            JLabel jLabel = new JLabel("Ground");
            jLabel.setBounds((min+max)/2 * 1200 / 6 + 100, 5 * 800 / 6, 60, 40);
            add(jLabel);
            super.paint(getGraphics());
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

        private void drawCCVS(CCVS ccvs) {
            int x1, y1, x2, y2;
            int node1Num = Integer.parseInt(ccvs.n1.name) - 1;
            int node2Num = Integer.parseInt(ccvs.n2.name) - 1;
            if (node1Num == -1)
                connectToEarth(node2Num, "VoltageCircleU.png", ccvs.NameH);
            else if (node2Num == -1)
                connectToEarth(node1Num, "VoltageCircleD.png", ccvs.NameH);
            else {
                x1 = node1Num % 6;
                y1 = 4 - node1Num / 6;
                x2 = node2Num % 6;
                y2 = 4 - node2Num / 6;
                if (x1 == x2) {
                    if (y1 > y2)
                        drawVertical(x1, y1, x2, y2, "VoltageCircleU.png", ccvs.NameH);
                    else
                        drawVertical(x1, y1, x2, y2, "VoltageCircleD.png", ccvs.NameH);
                } else if (y1 == y2) {
                    if (x1 > x2)
                        drawHorizontal(x1, y1, x2, y2, "VoltageCircleL.png", ccvs.NameH);
                    else
                        drawHorizontal(x1, y1, x2, y2, "VoltageCircleR.png", ccvs.NameH);
                }
            }
        }

        private void drawVCVS(VCVS vcvs) {
            int x1, y1, x2, y2;
            int node1Num = Integer.parseInt(vcvs.n1.name) - 1;
            int node2Num = Integer.parseInt(vcvs.n2.name) - 1;
            if (node1Num == -1)
                connectToEarth(node2Num, "VoltageCircleU.png", vcvs.NameE);
            else if (node2Num == -1)
                connectToEarth(node1Num, "VoltageCircleD.png", vcvs.NameE);
            else {
                x1 = node1Num % 6;
                y1 = 4 - node1Num / 6;
                x2 = node2Num % 6;
                y2 = 4 - node2Num / 6;
                if (x1 == x2) {
                    if (y1 > y2)
                        drawVertical(x1, y1, x2, y2, "VoltageCircleU.png", vcvs.NameE);
                    else
                        drawVertical(x1, y1, x2, y2, "VoltageCircleD.png", vcvs.NameE);
                } else if (y1 == y2) {
                    if (x1 > x2)
                        drawHorizontal(x1, y1, x2, y2, "VoltageCircleL.png", vcvs.NameE);
                    else
                        drawHorizontal(x1, y1, x2, y2, "VoltageCircleR.png", vcvs.NameE);
                }
            }
        }

        private void drawCCCS(CCCS cccs) {
            int x1, y1, x2, y2;
            int node1Num = Integer.parseInt(cccs.n1.name) - 1;
            int node2Num = Integer.parseInt(cccs.n2.name) - 1;
            if (node1Num == -1)
                connectToEarth(node2Num, "CurrentU.png", cccs.NameF);
            else if (node2Num == -1)
                connectToEarth(node1Num, "CurrentD.png", cccs.NameF);
            else {
                x1 = node1Num % 6;
                y1 = 4 - node1Num / 6;
                x2 = node2Num % 6;
                y2 = 4 - node2Num / 6;
                if (x1 == x2) {
                    if (y1 > y2)
                        drawVertical(x1, y1, x2, y2, "CurrentU.png", cccs.NameF);
                    else
                        drawVertical(x1, y1, x2, y2, "CurrentD.png", cccs.NameF);
                } else if (y1 == y2) {
                    if (x1 > x2)
                        drawHorizontal(x1, y1, x2, y2, "CurrentL.png", cccs.NameF);
                    else
                        drawHorizontal(x1, y1, x2, y2, "CurrentR.png", cccs.NameF);
                }
            }
        }

        private void drawVCCS(VCCS vccs) {
            int x1, y1, x2, y2;
            int node1Num = Integer.parseInt(vccs.n1.name) - 1;
            int node2Num = Integer.parseInt(vccs.n2.name) - 1;
            if (node1Num == -1)
                connectToEarth(node2Num, "CurrentU.png", vccs.NameG);
            else if (node2Num == -1)
                connectToEarth(node1Num, "CurrentD.png", vccs.NameG);
            else {
                x1 = node1Num % 6;
                y1 = 4 - node1Num / 6;
                x2 = node2Num % 6;
                y2 = 4 - node2Num / 6;
                if (x1 == x2) {
                    if (y1 > y2)
                        drawVertical(x1, y1, x2, y2, "CurrentU.png", vccs.NameG);
                    else
                        drawVertical(x1, y1, x2, y2, "CurrentD.png", vccs.NameG);
                } else if (y1 == y2) {
                    if (x1 > x2)
                        drawHorizontal(x1, y1, x2, y2, "CurrentL.png", vccs.NameG);
                    else
                        drawHorizontal(x1, y1, x2, y2, "CurrentR.png", vccs.NameG);
                }
            }
        }

        private void drawVS(voltageSource voltageSource1) {
            int x1, y1, x2, y2;
            int node1Num = Integer.parseInt(voltageSource1.n1.name) - 1;
            int node2Num = Integer.parseInt(voltageSource1.n2.name) - 1;
            if (voltageSource1.A != 0){
                if (node1Num == -1)
                    connectToEarth(node2Num, "AC.png", voltageSource1.NameV);
                else if (node2Num == -1)
                    connectToEarth(node1Num, "AC.png", voltageSource1.NameV);
                else {
                    x1 = node1Num % 6;
                    y1 = 4 - node1Num / 6;
                    x2 = node2Num % 6;
                    y2 = 4 - node2Num / 6;
                    if (x1 == x2) {
                        if (y1 > y2)
                            drawVertical(x1, y1, x2, y2, "AC.png", voltageSource1.NameV);
                        else
                            drawVertical(x1, y1, x2, y2, "AC.png", voltageSource1.NameV);
                    } else if (y1 == y2) {
                        if (x1 > x2)
                            drawHorizontal(x1, y1, x2, y2, "AC.png", voltageSource1.NameV);
                        else
                            drawHorizontal(x1, y1, x2, y2, "AC.png", voltageSource1.NameV);
                    }
                }
                return;
            }
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

        private void drawD(diode diode1) {
            int x1, y1, x2, y2;
            int node1Num = Integer.parseInt(diode1.n1.name) - 1;
            int node2Num = Integer.parseInt(diode1.n2.name) - 1;
            if (node1Num == -1)
                connectToEarth(node2Num, "DiodeU.png", diode1.NameD);
            else if (node2Num == -1)
                connectToEarth(node1Num, "DiodeD.png", diode1.NameD);
            else {
                x1 = node1Num % 6;
                y1 = 4 - node1Num / 6;
                x2 = node2Num % 6;
                y2 = 4 - node2Num / 6;
                if (x1 == x2) {
                    if (y1 > y2)
                        drawVertical(x1, y1, x2, y2, "DiodeU.png", diode1.NameD);
                    else
                        drawVertical(x1, y1, x2, y2, "DiodeD.png", diode1.NameD);
                } else if (y1 == y2) {
                    if (x1 > x2)
                        drawHorizontal(x1, y1, x2, y2, "DiodeL.png", diode1.NameD);
                    else
                        drawHorizontal(x1, y1, x2, y2, "DiodeR.png", diode1.NameD);
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
            int singleX = 0, singleY = 0;
            if (paralel2.get(string) == 1){
                singleY = 30;
            }
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
                jLabel.setBounds((int)((double)(x1 + x2) / 2.0 * (1200 / 6)) + 35, y1 * 800 / 6 + count * 40 + 12 + singleY, 70, 25);
            else if (address.contains("Current"))
                jLabel.setBounds((int)((double)(x1 + x2) / 2.0 * (1200 / 6)) + 35, y1 * 800 / 6 + count * 40 + 10 + singleY, 70, 32);
            else
                jLabel.setBounds((int)((double)(x1 + x2) / 2.0 * (1200 / 6)) + 35, y1 * 800 / 6 + count * 40 + singleY, 70, 50);
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
            jLabel1.setBounds((int)((double)(x1 + x2) / 2.0 * (1200 / 6)) + 60 + additionalX, y1 * 800 / 6 + count * 40 - 15 + additionalY + singleY, 70, 50);

            add(jLabel1);
            super.paint(getGraphics());

            Line2D line1 = new Line2D.Float(x1 * 1200 / 6 + 100, y1 * 800 / 6 + 80, (int)((double)(x1 + x2) / 2.0 * (1200 / 6)) + 45, y1 * 800 / 6 + count * 40 + 55 + singleY+8);
            Line2D line2 = new Line2D.Float((int)((double)(x1 + x2) / 2.0 * (1200 / 6)) + 110, y1 * 800 / 6 + count * 40 + 55 + singleY+8, x2 * 1200 / 6 + 100, y1 * 800 / 6 + 80);
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
            int singleX = 0, singleY = 0;
            if (paralel2.keySet().contains(string) && paralel2.get(string) == 1){
                singleX = 57;
            }
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
                jLabel.setBounds(x1 * 1200 / 6 + count * 40 + 22 + singleX, (int)((double)(y1 + y2) / 2.0 * (800 / 6)) - 5, 25, 45);
            else if (address.contains("Current"))
                jLabel.setBounds(x1 * 1200 / 6 + count * 40 + 15 + singleX, (int)((double)(y1 + y2) / 2.0 * (800 / 6)) - 5, 40, 50);
            else if (address.contains("Voltage"))
                jLabel.setBounds(x1 * 1200 / 6 + count * 40 + 17 + singleX, (int)((double)(y1 + y2) / 2.0 * (800 / 6)) - 5, 35, 50);
            else
                jLabel.setBounds(x1 * 1200 / 6 + count * 40 + singleX, (int)((double)(y1 + y2) / 2.0 * (800 / 6)) - 5, 70, 50);
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
            jLabel1.setBounds(x1 * 1200 / 6 + count * 40 + 8 + additionalX + singleX, (int)((double)(y1 + y2) / 2.0 * (800 / 6)) - 5 + additionalY, 70, 50);

            add(jLabel1);
            super.paint(getGraphics());
            Line2D line1 = new Line2D.Float(x1 * 1200 / 6 + 100, y1 * 800 / 6 + 80, x1 * 1200 / 6 + count * 40 + 42 + singleX, (int)((double)(y1 + y2) / 2.0 * (800 / 6)) + 27+8);
            lines.add(line1);
            Line2D line2 = new Line2D.Float(x1 * 1200 / 6 + count * 40 + 42 + singleX, (int)((double)(y1 + y2) / 2.0 * (800 / 6)) + 71+8, x2 * 1200 / 6 + 100, y2 * 800 / 6 + 80);
            lines.add(line2);

            paralel.put(string, count + 1);
        }

    }

    public static class chartPainterVoltage extends JFrame{
        branch x;
        double min,max;
        ArrayList<Double> VChart=new ArrayList<Double>(0);
        chartPainterVoltage(branch b) {
            x=b;
            for(int j=0;j<x.outputCurrent.size();j++){
                VChart.add(x.n1.outputVolt.get(j)-x.n2.outputVolt.get(j));
            }
            max=VChart.get(this.AndisMaxFinder(VChart));
            min=VChart.get(this.AndisMinFinder(VChart));
            setTitle("Voltage Chart");
            setSize(2000, 1000);
            setLayout(null);
            setVisible(true);
        }
        @Override
        public void paint(Graphics g){
            for(int j=0;j<VChart.size()-1;j++){
                g.drawLine(100+j*1700/VChart.size(),(int)(800-600*(VChart.get(j)-min)/(max-min)),100+(j+1)*1700/VChart.size(), (int)(800 -600 * (VChart.get(j + 1) - min) / (max - min)));

            }
            if(min<0&&max>0){
                g.drawLine(50,(int)(800+600*min/(max-min)),1850,(int)(800+600*min/(max-min)));
                g.drawLine(1830,(int)(800+600*min/(max-min))-20,1850,(int)(800+600*min/(max-min)));
                g.drawLine(1830,(int)(800+600*min/(max-min))+20,1850,(int)(800+600*min/(max-min)));
                g.drawLine(100,950,100,50);
                g.drawLine(80,70,100,50);
                g.drawLine(120,70,100,50);
                for (int j=1;j<10;j++){
                    g.drawString(Double.toString(min+(800-100*j)*(max-min)/600),50,100*j);
                }
                for (int j=1;j<19;j++){
                    g.drawString(Double.toString(100*(j-1)*VChart.size()/1700),100*j-20,900);
                }
            }
            else{
                g.drawLine(50,900,1850,900);
                g.drawLine(1830,880,1850,900);
                g.drawLine(1830,920,1850,900);
                g.drawLine(100,950,100,50);
                g.drawLine(80,70,100,50);
                g.drawLine(120,70,100,50);
                for (int j=1;j<10;j++){
                    g.drawString(Double.toString(min+(800-100*j)*(max-min)/600),50,100*j);
                }
                for (int j=1;j<19;j++){
                    g.drawString(Double.toString(100*(j-1)*VChart.size()/1700),100*j-20,900);
                }
            }

        }
        public int AndisMinFinder(ArrayList<Double> A){
            int i=0;
            for(int j=0;j<A.size();j++)
                if(A.get(j)<A.get(i))
                    i=j;
            return i;
        }
        public int AndisMaxFinder(ArrayList<Double> A){
            int i=0;
            for(int j=0;j<A.size();j++)
                if(A.get(j)>A.get(i))
                    i=j;
            return i;
        }
    }
    public static class chartPainterCurrent extends JFrame{
        branch x;
        double max,min;
        chartPainterCurrent(branch b) {
            x=b;
            max=x.outputCurrent.get(this.AndisMaxFinder(x.outputCurrent));
            min=x.outputCurrent.get(this.AndisMinFinder(x.outputCurrent));
            setTitle("Current Chart");
            setSize(2000, 1000);
            setLayout(null);
            setVisible(true);
        }
        @Override
        public void paint(Graphics g){
            for(int j=0;j<x.outputCurrent.size()-1;j++){
                g.drawLine(100+j*1700/x.outputCurrent.size(),(int)(800-600*(x.outputCurrent.get(j)-min)/(max-min)),100+(j+1)*1700/x.outputCurrent.size(), (int)(800 -600 * (x.outputCurrent.get(j + 1) - min) / (max - min)));

            }
            if(min<0&&max>0){
                g.drawLine(50,(int)(800+600*min/(max-min)),1850,(int)(800+600*min/(max-min)));
                g.drawLine(1830,(int)(800+600*min/(max-min))-20,1850,(int)(800+600*min/(max-min)));
                g.drawLine(1830,(int)(800+600*min/(max-min))+20,1850,(int)(800+600*min/(max-min)));
                g.drawLine(100,950,100,50);
                g.drawLine(80,70,100,50);
                g.drawLine(120,70,100,50);
                for (int j=1;j<10;j++){
                    g.drawString(Double.toString(min+(800-100*j)*(max-min)/600),50,100*j);
                }
                for (int j=1;j<19;j++){
                    g.drawString(Double.toString(100*(j-1)*x.outputCurrent.size()/1700),100*j-20,900);
                }
            }
            else{
                g.drawLine(50,900,1850,900);
                g.drawLine(1830,880,1850,900);
                g.drawLine(1830,920,1850,900);
                g.drawLine(100,950,100,50);
                g.drawLine(80,70,100,50);
                g.drawLine(120,70,100,50);
                for (int j=1;j<10;j++){
                    g.drawString(Double.toString(min+(800-100*j)*(max-min)/600),50,100*j);
                }
                for (int j=1;j<19;j++){
                    g.drawString(Double.toString(100*(j-1)*x.outputCurrent.size()/1700),100*j-20,900);
                }
            }

        }
        public int AndisMinFinder(ArrayList<Double> A){
            int i=0;
            for(int j=0;j<A.size();j++)
                if(A.get(j)<A.get(i))
                    i=j;
            return i;
        }
        public int AndisMaxFinder(ArrayList<Double> A){
            int i=0;
            for(int j=0;j<A.size();j++)
                if(A.get(j)>A.get(i))
                    i=j;
            return i;
        }
    }
    public static class chartPainterPower extends JFrame{
        branch x;
        double min,max;
        ArrayList<Double> PChart=new ArrayList<Double>(0);
        chartPainterPower(branch b) {
            x=b;
            for(int j=0;j<x.outputCurrent.size();j++){
                PChart.add((x.n1.outputVolt.get(j)-x.n2.outputVolt.get(j))*x.outputCurrent.get(j));
            }
            max=PChart.get(this.AndisMaxFinder(PChart));
            min=PChart.get(this.AndisMinFinder(PChart));
            setTitle("Power Chart");
            setSize(2000, 1000);
            setLayout(null);
            setVisible(true);
        }
        @Override
        public void paint(Graphics g){
            for(int j=0;j<PChart.size()-1;j++){
                g.drawLine(100+j*1700/PChart.size(),(int)(800-600*(PChart.get(j)-min)/(max-min)),100+(j+1)*1700/PChart.size(), (int)(800 -600 * (PChart.get(j + 1) - min) / (max - min)));

            }
            if(min<0&&max>0){
                g.drawLine(50,(int)(800+600*min/(max-min)),1850,(int)(800+600*min/(max-min)));
                g.drawLine(1830,(int)(800+600*min/(max-min))-20,1850,(int)(800+600*min/(max-min)));
                g.drawLine(1830,(int)(800+600*min/(max-min))+20,1850,(int)(800+600*min/(max-min)));
                g.drawLine(100,950,100,50);
                g.drawLine(80,70,100,50);
                g.drawLine(120,70,100,50);
                for (int j=1;j<10;j++){
                    g.drawString(Double.toString(min+(800-100*j)*(max-min)/600),50,100*j);
                }
                for (int j=1;j<19;j++){
                    g.drawString(Double.toString(100*(j-1)*PChart.size()/1700),100*j-20,900);
                }
            }
            else{
                g.drawLine(50,900,1850,900);
                g.drawLine(1830,880,1850,900);
                g.drawLine(1830,920,1850,900);
                g.drawLine(100,950,100,50);
                g.drawLine(80,70,100,50);
                g.drawLine(120,70,100,50);
                for (int j=1;j<10;j++){
                    g.drawString(Double.toString(min+(800-100*j)*(max-min)/600),50,100*j);
                }
                for (int j=1;j<19;j++){
                    g.drawString(Double.toString(100*(j-1)*PChart.size()/1700),100*j-20,900);
                }
            }

        }

        public int AndisMinFinder(ArrayList<Double> A){
            int i=0;
            for(int j=0;j<A.size();j++)
                if(A.get(j)<A.get(i))
                    i=j;
            return i;
        }
        public int AndisMaxFinder(ArrayList<Double> A){
            int i=0;
            for(int j=0;j<A.size();j++)
                if(A.get(j)>A.get(i))
                    i=j;
            return i;
        }
    }


    public static class ConsolMatni extends JFrame{
        String x;
        double[] T=new double[2];
        {T[0]=-1;T[1]=-1;}
        ConsolMatni() throws IOException {
            JButton Load=new JButton("Load");
            JButton Run=new JButton("Run");
            JTextArea Consol=new JTextArea();
            Consol.setEditable(true);
            Consol.setLineWrap(true);
            Consol.setWrapStyleWord(true);
            JLabel CONSOLMATNI=new JLabel("consol matni :");
            Load.setBounds(0,0,200,50);
            Run.setBounds(200,0,200,50);
            Consol.setBounds(200,100,600,300);
            CONSOLMATNI.setBounds(100,100,100,50);



            Scanner scanner=new Scanner(file);
            while (scanner.hasNextLine()){

                String line;
                line=scanner.nextLine();
                Consol.append(line+"\n");
            }

            final FileWriter[] fileWriter = {new FileWriter(file, false)};


            Load.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    Consol.setText("");
                    JFrame f1 = new JFrame("Testing File Selection");
                    FileSystemView fsv;
                    fsv = FileSystemView.getFileSystemView();
                    file = new File("C:\\");
                    JFileChooser fileChooser = new JFileChooser(file,fsv);
                    int response = fileChooser.showOpenDialog(f1);
                    if (response == JFileChooser.APPROVE_OPTION){
                        file=fileChooser.getSelectedFile();
                        try {
                            Scanner scanner = new Scanner(file);
                            //System.out.println("asda");
                            while (scanner.hasNextLine()) {

                                String line;
                                line = scanner.nextLine();
                                Consol.append(line + "\n");
                            }
                        }
                        catch (Exception z){

                        }
                    }
                }
                });


            Run.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    try{
                        x=Consol.getText();
                        fileWriter[0] =new FileWriter(file,false);
                        fileWriter[0].write(x);
                        fileWriter[0].close();
                    }
                    catch (Exception w){
                        System.out.println(w);
                    }
                    input=new ArrayList<String>(0);
                    N=new ArrayList<node>(0);
                    U=new ArrayList<nodes>(0);
                    R=new ArrayList<resistor>(0);
                    CS=new ArrayList<currentSource>(0);
                    VS=new ArrayList<voltageSource>(0);
                    C=new ArrayList<capacitor>(0);
                    L=new ArrayList<inductor>(0);
                    G=new ArrayList<VCCS>(0);
                    F=new ArrayList<CCCS>(0);
                    E=new ArrayList<VCVS>(0);
                    H=new ArrayList<CCVS>(0);
                    D=new ArrayList<diode>(0);


                    T=tahlilMadar(file);


                    int error=1;
                    Scanner scanner =new Scanner(System.in);
                    String s=scanner.nextLine();
                    s = s.trim();
                    s = s.replaceAll("( )+", " ");
                    while (!s.equals("END")) {
                        try {
                            error = consoleInput(s, T[1], T[0]);
                            s += 1 / error;
                            s = scanner.nextLine();
                            s = s.trim();
                            s = s.replaceAll("( )+", " ");
                            if(s.equals("END")){
                                break;
                            }
                        }
                        catch (Exception q){
                            s="";
                            System.out.println("ERROR");
                            s=scanner.nextLine();
                            s = s.trim();
                            s = s.replaceAll("( )+", " ");
                            if(s.equals("END")){
                                break;
                            }
                        }
                    }
                    System.out.println("programm ends");
                }
            });

            add(Load);
            add(Run);
            add(Consol);
            add(CONSOLMATNI);




            setTitle("CONSOL MATNI");
            setSize(1000, 500);
            setLayout(null);
            setVisible(true);
        }
    }


    public static double[] tahlilMadar(File file){
        node n;
        resistor r;
        currentSource cs;
        voltageSource vs;
        capacitor c;
        inductor l;
        diode d;
        CCCS cccs;
        VCCS vccs;
        VCVS vcvs;
        CCVS ccvs;
        double T=-1, dT=-1,dV=-1, dI=-1;

        double[] tout=new double[2];
        {tout[0]=-1;tout[1]=-1;}

        int lineNumber=1, errorType=0;
        try {

            File output = new File("output.txt");
            FileWriter fileWriter = new FileWriter(output);
            Scanner sc = new Scanner(file);

            String s = sc.nextLine();
            s = s.trim();
            s = s.replaceAll("( )+", " ");


            while (s.indexOf(".tran") == -1) {
                if (s.charAt(0) == '*') {}
                else if (s.charAt(0) == 'R') {
                    r = new resistor(s);
                    r.outputCurrent.add(0.0);
                    R.add(r);
                    input.add("R");
                    if (r.R < 0) {
                        s += 0 / 0;
                    }
                }
                else if (s.charAt(0) == 'L') {
                    l = new inductor(s);
                    l.outputCurrent.add(0.0);
                    L.add(l);
                    input.add("L");
                    if (l.L < 0) {
                        s += 0 / 0;
                    }
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
                    if (c.C < 0) {
                        s += 0 / 0;
                    }
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
                else if (s.charAt(0) == 'E') {
                    vcvs = new VCVS(s);
                    E.add(vcvs);
                    input.add("E");
                }
                else if (s.charAt(0) == 'H') {
                    ccvs = new CCVS(s);
                    H.add(ccvs);
                    input.add("H");
                }
                else if (s.charAt(0) == 'D') {
                    d = new diode(s);
                    D.add(d);
                    input.add("D");
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
                    if (dT <= 0) {
                        s += 0 / 0;
                    }
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
                    if (dV <= 0) {
                        s += 0 / 0;
                    }
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
                    if (dI <= 0) {
                        s += 0 / 0;
                    }
                }
                else {
                    s += 0 / 0;
                }
                lineNumber++;
                s = sc.nextLine();
                s = s.trim();
                s = s.replaceAll("( )+", " ");
            }

            s = s.substring(s.indexOf(" ") + 1);
            s = s.replaceAll("k", "000");
            s = s.replaceAll("M", "000000");
            s = s.replaceAll("G", "000000000");
            if (s.indexOf("m") != -1) {
                T = 0.001;
                s = s.replaceAll("m", "");
            }
            else if (s.indexOf("u") != -1) {
                T = 0.000001;
                s = s.replaceAll("u", "");
            }
            else if (s.indexOf("n") != -1) {
                T = 0.000000001;
                s = s.replaceAll("n", "");
            }
            else if (s.indexOf("p") != -1) {
                T = 0.000000000001;
                s = s.replaceAll("p", "");
            }
            else {
                T = 1;
            }
            T *= Double.parseDouble(s);
            if (T <= 0) {
                s += 0 / 0;
            }



            for (int j = 0; j < U.size(); j++) {
                U.get(j).n.get(0).outputVolt.add(0.0);
            }
            for (int j = 0; j < U.size(); j++) {
                U.get(j).addInputUnionVolts0(0, 0, 0, dT);
                /*for (int x = 0; x < U.get(j).n.size(); x++) {
                    System.out.println(U.get(j).n.get(x).name + " : " + U.get(j).n.get(x).union + " in: " + U.get(j).union + " Volt: " + U.get(j).n.get(x).outputVolt);
                }*/
            }

            int checkingGround = searchNode("0");
            if (checkingGround == -1) {
                errorType = -4;
                s += 0 / 0;
            }
            int checkingVS = 0;

            if (dT > 0 && T > 0 && dI > 0 && dV > 0) {
                for (int i = 0; i <= T / dT; i += 1) {
                    checkingGround = checkGround(dT, i);
                    if (checkingGround == -1) {
                        errorType = -4;
                        s += 0 / 0;
                    }
                    checkingVS = checkVS(dT, i);
                    if (checkingVS == -1) {
                        errorType = -3;
                        s += 0 / 0;
                    }
                    if (checkCS(dT, i) == -1) {
                        errorType = -2;
                        s += 0 / 0;
                    }
                    calcNodeVolts(dT, dV, dI, i);
                    calcBranchCurrents(dT, dV, dI, i);
                    /*System.out.println("RAN : "+i);
                    for (int j = 0; j < U.size(); j++) {
                        for (int x = 0; x < U.get(j).n.size(); x++) {
                            System.out.println(U.get(j).n.get(x).name + " : " + U.get(j).n.get(x).union + " in: " + U.get(j).union + " Volt: " + U.get(j).n.get(x).outputVolt);
                        }
                    }*/
                }
            }
            else {
                errorType = -1;
                s += 0 / 0;
            }


            chapOutput(fileWriter);
            fileWriter.close();


            graphProject GraphProject = new graphProject();

            tout[0] = T;
            tout[1] = dT;

        }
        catch (Exception e){
            System.out.println(e);
            if(errorType==0) {
                System.out.println("Error on line: " + lineNumber);
            }
            else{
                System.out.println("Error : "+errorType);
            }
        }
        return tout;
    }


    public static void main(String[] args) {
        try {

            new ConsolMatni();



        }
        catch (Exception e){
            System.out.println(e);
        }
    }
}
