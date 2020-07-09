import javax.swing.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class project {

    static ArrayList<node> N=new ArrayList<node>(0);
    static ArrayList<resistor> R=new ArrayList<resistor>(0);
    static ArrayList<currentSource> CS=new ArrayList<currentSource>(0);
    static ArrayList<voltageSource> VS=new ArrayList<voltageSource>(0);
    static ArrayList<capacitor> C=new ArrayList<capacitor>(0);
    static ArrayList<inductor> L=new ArrayList<inductor>(0);
    static ArrayList<currentControledCurrentSource> CCCS=new ArrayList<currentControledCurrentSource>(0);
    static ArrayList<removedNode> RemovedNode=new ArrayList<removedNode>(0);



    public static int searchNode(String s){
        for(int i=0;i<N.size();i++){
            if(N.get(i).name.equals(s)){
                return i;
            }
        }
        return -1;
    }
    public static int searchRemovedNode(String s){
        for(int i=0;i<RemovedNode.size();i++){
            if(RemovedNode.get(i).name.equals(s)){
                return i;
            }
        }
        return -1;
    }


    public static class node{
        String name, output;
        double volt;
        node(String s){
            name=s;
        }
    }
    public static class removedNode extends node{
        String n1name;
        removedNode(String s) {
            super(s);
        }
    }
    abstract public static class branch{
        String output;
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
        void calc_I(){
            I=(n1.volt-n2.volt)/R;
        }
        void calc_v1(int t){
            n1.volt=n2.volt+I*R;
        }
        void calc_v2(int t){
            n2.volt=n1.volt-I*R;
        }
    }
    public static class capacitor extends branch {
        String NameC;
        double C;
        capacitor(String s){
            NameC=s.substring(0,s.indexOf(" "));
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
    public static class currentControledCurrentSource extends branch {
        String NameI, NameEleman;
        double a;
        currentControledCurrentSource(String s) {
            NameI = s.substring(0, s.indexOf(" "));
            s = s.substring(s.indexOf(" ") + 1);
            node n = new node(s.substring(0, s.indexOf(" ")));
            int i = searchNode(s.substring(0, s.indexOf(" ")));
            if (i == -1) {
                N.add(n);
                n1 = n;
            }
            else {
                n1 = N.get(i);
            }
            s = s.substring(s.indexOf(" ") + 1);
            n = new node(s.substring(0, s.indexOf(" ")));
            i = searchNode(s.substring(0, s.indexOf(" ")));
            if (i == -1) {
                N.add(n);
                n2 = n;
            }
            else {
                n2 = N.get(i);
            }
            s = s.substring(s.indexOf(" ") + 1);
            NameEleman=s.substring(0, s.indexOf(" "));
            s=s.substring(s.indexOf(" ")+1);
            a=Double.parseDouble(s);
        }
        public void update(){
            if (NameEleman.contains("C")) {
                for (int j = 0; j < C.size(); j++) {
                    if (C.get(j).NameC.equals(NameEleman)) {
                        I = a * C.get(j).I;
                    }
                }
            }
            else if (NameEleman.contains("R")) {
                for (int j = 0; j < R.size(); j++) {
                    if (R.get(j).NameR.equals(NameEleman)) {
                        I = a* R.get(j).I;
                    }
                }
            }
            else if (NameEleman.contains("I")) {
                for (int j = 0; j < CS.size(); j++) {
                    if (CS.get(j).NameI.equals(NameEleman)) {
                        I =a*CS.get(j).I;
                    }
                }
            }
            else if (NameEleman.contains("L")) {
                for (int j = 0; j < L.size(); j++){
                    if (L.get(j).NameL.equals(NameEleman)) {
                        I = a* L.get(j).I;
                    }
                }
        }
        }
    }

    static class Gauss_Jordan_Elimination {
        private static final double EPSILON = 1e-8;
        private static int N;
        private static double[][] a;
        public Gauss_Jordan_Elimination(double[][] A, double[] b) {
            N = b.length;
            a = new double[N][N+N+1];
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    a[i][j] = A[i][j];
            for (int i = 0; i < N; i++)
                a[i][N+i] = 1.0;
            for (int i = 0; i < N; i++)
                a[i][N+N] = b[i];
            show();
            solve();
            assert check(A, b);
        }
        private void solve() {
            for (int p = 0; p < N; p++) {
                int max = p;
                for (int i = p+1; i < N; i++) {
                    if (Math.abs(a[i][p]) > Math.abs(a[max][p])) {
                        max = i;
                    }
                }
                swap(p, max);
                if (Math.abs(a[p][p]) <= EPSILON) {
                    continue;
                }
                pivot(p, p);
            }
        }
        private void swap(int row1, int row2) {
            double[] temp = a[row1];
            a[row1] = a[row2];
            a[row2] = temp;
        }
        private void pivot(int p, int q) {
            for (int i = 0; i < N; i++) {
                double alpha = a[i][q] / a[p][q];
                for (int j = 0; j <= N+N; j++) {
                    if (i != p && j != q) a[i][j] -= alpha * a[p][j];
                }
            }
            for (int i = 0; i < N; i++)
                if (i != p) a[i][q] = 0.0;
            for (int j = 0; j <= N+N; j++)
                if (j != q) a[p][j] /= a[p][q];
            a[p][q] = 1.0;
        }
        public double[] primal() {
            double[] x = new double[N];
            for (int i = 0; i < N; i++) {
                if (Math.abs(a[i][i]) > EPSILON)
                    x[i] = a[i][N+N] / a[i][i];
                else if (Math.abs(a[i][N+N]) > EPSILON)
                    return null;
            }
            return x;
        }
        public double[] dual() {
            double[] y = new double[N];
            for (int i = 0; i < N; i++) {
                if ( (Math.abs(a[i][i]) <= EPSILON) && (Math.abs(a[i][N+N]) > EPSILON) ) {
                    for (int j = 0; j < N; j++)
                        y[j] = a[i][N+j];
                    return y;
                }
            }
            return null;
        }
        public boolean isFeasible() {
            return primal() != null;
        }
        private static void show() {
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    System.out.print(a[i][j]+"  ");
                }
                System.out.print("|");
                for (int j = N; j < N+N; j++) {
                    System.out.print(a[i][j]+"  ");
                }
                System.out.println("| "+a[i][N+N]);
            }
            System.out.println();
        }
        private boolean check(double[][] A, double[] b) {
            if (isFeasible()) {
                double[] x = primal();
                for (int i = 0; i < N; i++) {
                    double sum = 0.0;
                    for (int j = 0; j < N; j++) {
                        sum += A[i][j] * x[j];
                    }
                    if (Math.abs(sum - b[i]) > EPSILON) {
                        System.out.println("not feasible");
                        System.out.println(i+" = "+b[i]+", sum = "+sum+"\n");
                        return false;
                    }
                }
                return true;
            }
            else {
                double[] y = dual();
                for (int j = 0; j < N; j++) {
                    double sum = 0.0;
                    for (int i = 0; i < N; i++) {
                        sum += A[i][j] * y[i];
                    }
                    if (Math.abs(sum) > EPSILON) {
                        System.out.println("invalid certificate of infeasibility");
                        System.out.println("sum = "+sum+"\n");
                        return false;
                    }
                }
                double sum = 0.0;
                for (int i = 0; i < N; i++) {
                    sum += y[i] * b[i];
                }
                if (Math.abs(sum) < EPSILON) {
                    System.out.println("invalid certificate of infeasibility");
                    System.out.println("yb  = "+sum+"\n");
                    return false;
                }
                return true;
            }
        }
        public static void test(double[][] A, double[] b) {
            Gauss_Jordan_Elimination gaussian = new Gauss_Jordan_Elimination(A, b);
            if (gaussian.isFeasible()) {
                System.out.println("Solution to Ax = b");
                double[] x = gaussian.primal();
                for (int i = 0; i < x.length; i++) {
                    System.out.println(" "+x[i]+"\n");
                }
            }
            else {
                System.out.println("Certificate of infeasibility");

                double[] y = gaussian.dual();
                for (int j = 0; j < y.length; j++) {
                    System.out.println(" "+y[j]+"\n");
                }

            }
            updateNodeVolts();
            show();
            System.out.println();
        }
        public static void updateNodeVolts(){
            for(int i=0;i<a.length;i++){
                project.N.get(i).volt=a[i][a[i].length-1];
            }
        }
    }


    static double calcGii(int i){
        double gii=0;
        for(int k=0;k<R.size();k++){
            if(R.get(k).n2.name.equals(N.get(i).name)||R.get(k).n1.name.equals(N.get(i).name)){
                if(!R.get(k).n2.name.equals(R.get(k).n1.name))
                    gii+=1.0/(R.get(k).R);
            }
        }
        return gii;
    }
    static double calcGij(int i, int j){
        double gij=0;
        for(int k=0;k<R.size();k++){
            if(R.get(k).n2.name.equals(N.get(i).name)&&R.get(k).n1.name.equals(N.get(j).name)){
                gij-=1.0/(R.get(k).R);
            }
            else if(R.get(k).n1.name.equals(N.get(i).name)&&R.get(k).n2.name.equals(N.get(j).name)){
                gij-=1.0/(R.get(k).R);
            }
        }
        return gij;
    }
    static double calcJi(int i){
        double ji=0;
        for(int k=0;k<CS.size();k++){
            if(CS.get(k).n2.name.equals(N.get(i).name)){
                ji+=CS.get(k).I;
            }
            if(CS.get(k).n1.name.equals(N.get(i).name)){
                ji-=CS.get(k).I;
            }
        }
        return ji;
    }




    public static void updateMadar(double deltat, double t){
        Pattern pattern=Pattern.compile("(.+)_(.+)");
        Pattern pattern1=Pattern.compile("IL(.+)");
        Matcher matcher1;
        Matcher matcher2;
        double I0=0;
        int a=0;
        for(int i=0;i<C.size();i++) {
            a=0;
            I0=0;
            for(int j=0;j<CS.size();j++) {
                matcher1=pattern.matcher(CS.get(j).NameI);
                matcher1.find();
                if(CS.get(j).NameI.contains("IVC"+Integer.toString(i))){
                    if (a==0) {
                        for (int k = 0; k < CS.size(); k++) {
                            if (CS.get(k).NameI.contains("IVC" + Integer.toString(i))) {
                                I0 = CS.get(k).I;
                                if(R.get(Integer.parseInt(matcher1.group(2))).n1.name.equals(CS.get(k).n1.name)){
                                    I0+=(R.get(Integer.parseInt(matcher1.group(2))).n1.volt-R.get(Integer.parseInt(matcher1.group(2))).n2.volt)/R.get(Integer.parseInt(matcher1.group(2))).R;
                                }
                                else if(R.get(Integer.parseInt(matcher1.group(2))).n2.name.equals(CS.get(k).n1.name)) {
                                    I0+=(R.get(Integer.parseInt(matcher1.group(2))).n2.volt - R.get(Integer.parseInt(matcher1.group(2))).n1.volt) / R.get(Integer.parseInt(matcher1.group(2))).R;
                                }
                                a = 1;
                            }
                        }
                    }
                    C.get(i).n2.volt+=I0*deltat/C.get(i).C;
                    CS.get(j).I+=((deltat/C.get(i).C)/R.get(Integer.parseInt(matcher1.group(2))).R)*I0;
                    System.out.println("Cvoltage : "+C.get(i).n2.name+" : "+C.get(i).n2.volt);
                }
            }
            C.get(i).I=I0;
        }
        for(int i=0;i<CS.size();i++){
            if(CS.get(i).NameI.contains("IL")) {
                matcher2=pattern1.matcher(CS.get(i).NameI);
                matcher2.find();
                System.out.println(L.get(Integer.parseInt(matcher2.group(1))).n1.name+"->"+L.get(Integer.parseInt(matcher2.group(1))).n2.name);
                CS.get(i).I+=(deltat/L.get(Integer.parseInt(matcher2.group(1))).L)*(CS.get(i).n1.volt-CS.get(i).n2.volt);
                L.get(Integer.parseInt(matcher2.group(1))).I=CS.get(i).I;

            }
            else if(CS.get(i).A!=0){
                CS.get(i).I += CS.get(i).A * (Math.sin(CS.get(i).w * ((t+1) * deltat) + CS.get(i).p)) - CS.get(i).A * (Math.sin(CS.get(i).w * (t*deltat) + CS.get(i).p));
                System.out.println("Iin: : "+CS.get(i).A * (Math.sin(CS.get(i).w * ((t+1) * deltat))));
                System.out.println("Iin: : "+CS.get(i).A * (Math.sin(CS.get(i).w * (t*deltat))));
            }
            else if(CS.get(i).NameI.contains("IV")){
                for (int j=0, k=searchNode(VS.get(j).n2.name);j<VS.size();j++){
                    if(CS.get(i).NameI.contains("I"+VS.get(j).NameV)){
                        CS.get(i).I/=(VS.get(j).V);
                        if(VS.get(j).A!=0) {
                            RemovedNode.get(k).volt += VS.get(j).A * Math.sin(VS.get(j).w * ((t + 1) * deltat) + VS.get(j).p) - VS.get(j).A * Math.sin(VS.get(j).w * (t * deltat) + VS.get(j).p);
                        }
                        VS.get(j).V += VS.get(j).A * Math.sin(VS.get(j).w * ((t+1) * deltat) + VS.get(j).p) - VS.get(j).A * Math.sin(VS.get(j).w * (t * deltat) + VS.get(j).p);
                        CS.get(i).I*=(VS.get(j).V);
                    }
                }
            }
        }
        for(int i=0;i<R.size();i++){
            R.get(i).calc_I();
        }
        for(int i=0;i<CCCS.size();i++){
            CCCS.get(i).update();
        }
        for(int k=0;k<VS.size();k++){
            if(VS.get(k).NameV.contains("C")){
                for(int g=0;g<RemovedNode.size();g++){
                    if (VS.get(k).n2.name.equals(RemovedNode.get(g).name)){
                        for(int m=0;m<C.size();m++){
                            if(VS.get(k).NameV.contains(C.get(m).NameC))
                                RemovedNode.get(g).volt=C.get(m).n2.volt+C.get(m).n1.volt;
                        }
                    }
                }
            }
            else{
                for (int g=0;g<RemovedNode.size();g++){
                    if (VS.get(k).n2.name.equals(RemovedNode.get(g).name)){
                        RemovedNode.get(g).volt=VS.get(k).n1.volt+VS.get(k).V;
                    }
                }
            }
        }
    }


    public static void replaceL(int i){
        String n2VName=new String(L.get(i).n2.name);
        String n1VName=new String(L.get(i).n1.name);
        int j1=searchNode(n1VName), j2=searchNode(n2VName);
        currentSource cs=new currentSource("IL"+Integer.toString(i)+" "+n1VName+" "+n2VName+" "+ Double.toString(L.get(i).I)+" 0 0 0");
        CS.add(cs);
    }
    public static void replaceVS(int i){
        String n2VName=new String(VS.get(i).n2.name);
        String n1VName=new String(VS.get(i).n1.name);
        int j=searchNode(n2VName);
        for(int k=0;k<R.size();k++){
            if(R.get(k).n1.name.equals(n2VName)){
                System.out.println(n2VName);
                currentSource csV=new currentSource("I"+VS.get(i).NameV+"_"+Integer.toString(k)+" "+n1VName+" "+R.get(k).n2.name+" "+Double.toString(VS.get(i).V/R.get(k).R)+" 0 0 0");
                System.out.println("1 "+R.get(k).NameR+" I"+VS.get(i).NameV+"_"+Integer.toString(k)+" "+n1VName+" "+R.get(k).n2.name+" "+Double.toString(VS.get(i).V/R.get(k).R)+" 0 0 0");
                CS.add(csV);
            }
            else if(R.get(k).n2.name.equals(n2VName)){
                System.out.println(n2VName);
                currentSource csV=new currentSource("I"+VS.get(i).NameV+"_"+Integer.toString(k)+" "+n1VName+" "+R.get(k).n1.name+" "+Double.toString(VS.get(i).V/R.get(k).R)+" 0 0 0");
                System.out.println("2 "+R.get(k).NameR+" I"+VS.get(i).NameV+"_"+Integer.toString(k)+" "+n1VName+" "+R.get(k).n1.name+" "+Double.toString(VS.get(i).V/R.get(k).R)+" 0 0 0");
                CS.add(csV);
            }
        }
        for (int k=0;k<R.size();k++){
            if(R.get(k).n1.name.equals(n2VName)){
                R.get(k).n1.name=n1VName;
            }
            else if(R.get(k).n2.name.equals(n2VName)){
                R.get(k).n2.name=n1VName;
            }
        }
        removedNode rn=new removedNode(n2VName);
        rn.n1name=n1VName;
        rn.volt=VS.get(i).V+N.get(searchNode(n1VName)).volt;
        RemovedNode.add(rn);
        N.remove(j);
    }
    public static void replaceC(int i){
        String n2VName=new String(C.get(i).n2.name);
        String n1VName=new String(C.get(i).n1.name);
        int j1=searchNode(n1VName), j2=searchNode(n2VName);
        voltageSource vc=new voltageSource("VC"+Integer.toString(i)+" "+n1VName+" "+n2VName+" "+ Double.toString(C.get(i).n2.volt-C.get(i).n1.volt)+" 0 0 0");
        VS.add(vc);
    }



    public static void main(String[] args) {
        resistor r;
        currentSource cs;
        voltageSource vs;
        capacitor c;
        inductor l;
        currentControledCurrentSource cccs;
        double T=1, dT=1;
        Scanner sc=new Scanner(System.in);
        String s=sc.nextLine();
        s=s.trim();
        s=s.replaceAll("( )+", " ");

        while (!s.equals(".end")){
            if(s.charAt(0)=='R'){
                r=new resistor(s);
                R.add(r);
            }
            else if(s.charAt(0)=='L'){
                l=new inductor(s);
                L.add(l);
            }
            else if(s.charAt(0)=='I'){
                cs=new currentSource(s);
                CS.add(cs);
            }
            else if(s.charAt(0)=='V'){
                vs=new voltageSource(s);
                VS.add(vs);
            }
            else if(s.charAt(0)=='C') {
                c=new capacitor(s);
                C.add(c);
            }
            else if(s.charAt(0)=='F'){
                cccs=new currentControledCurrentSource(s);
                CCCS.add(cccs);
            }
            else if(s.indexOf(".tran")!=-1){
                s=s.substring(s.indexOf(" ")+1);
                s=s.replaceAll("k", "000");
                s=s.replaceAll("M", "000000");
                s=s.replaceAll("G", "000000000");
                if(s.indexOf("m")!=-1){
                    T=0.001;
                    s=s.replaceAll("m","");
                }
                else if(s.indexOf("u")!=-1){
                    T=0.000001;
                    s=s.replaceAll("u", "");
                }
                else if(s.indexOf("n")!=-1){
                    T=0.000000001;
                    s=s.replaceAll("n", "");
                }
                else if(s.indexOf("p")!=-1){
                    T=0.000000000001;
                    s=s.replaceAll("p", "");
                }
                else {
                    T=1;
                }
                T*=Double.parseDouble(s);
            }
            else if(s.indexOf("dT")!=-1){
                s=s.substring(s.indexOf(" ")+1);
                s=s.replaceAll("k", "000");
                s=s.replaceAll("M", "000000");
                s=s.replaceAll("G", "000000000");
                if(s.indexOf("m")!=-1){
                    dT=0.001;
                    s=s.replaceAll("m","");
                }
                else if(s.indexOf("u")!=-1){
                    dT=0.000001;
                    s=s.replaceAll("u", "");
                }
                else if(s.indexOf("n")!=-1){
                    dT=0.000000001;
                    s=s.replaceAll("n", "");
                }
                else if(s.indexOf("p")!=-1){
                    dT=0.000000000001;
                    s=s.replaceAll("p", "");
                }
                else {
                    dT=1;
                }
                dT*=Double.parseDouble(s);

            }


            s=sc.nextLine();
            s=s.trim();
            s=s.replaceAll("( )+", " ");
        }


        for(int i=0;i<L.size();i++) {
            replaceL(i);
        }
        for(int i=0;i<C.size();i++) {
            replaceC(i);
        }
        for(int i=VS.size()-1;i>-1;i--) {
            replaceVS(i);
        }

        for(int t=0;t<T/dT;t++) {
            if (CS.size() > 0) {
                int n = N.size() - 1;
                double[][] mat = new double[n][n];
                double[] constants = new double[n];
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        if (i == j) {
                            mat[i][j] = calcGii(i);
                        }
                        else {
                            mat[i][j] = calcGij(i, j);
                        }
                    }
                    constants[i] = calcJi(i);
                }
                Gauss_Jordan_Elimination.test(mat, constants);
                for(int x=0;x<N.size();x++){
                    System.out.println(N.get(x).name+"  "+N.get(x).volt);
                }
                for (int x=0;x<RemovedNode.size();x++){
                    System.out.println(RemovedNode.get(x).name+"  "+Double.toString(RemovedNode.get(x).volt+N.get(searchNode(RemovedNode.get(x).n1name)).volt));
                }
            }

            if (dT != 1) {
                updateMadar(dT, t);
            }
        }
    }
}
