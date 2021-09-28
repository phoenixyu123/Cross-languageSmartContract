import it.unisa.dia.gas.jpbc.Element;//注意这个Java之中有一个同名的接口Element，一定要选这个
import it.unisa.dia.gas.jpbc.Field;
import it.unisa.dia.gas.jpbc.Pairing;
import it.unisa.dia.gas.jpbc.PairingParameters;
import it.unisa.dia.gas.plaf.jpbc.pairing.PairingFactory;
import it.unisa.dia.gas.plaf.jpbc.pairing.a.TypeACurveGenerator;

import java.math.BigInteger;
import java.util.Scanner;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * 提供随机输入
 * 【注意！】涉及双线性对运算的（pairing），建议线程不超过8个
 * @Author 北方工业大学 Yigang Yang 1565636387@qq.com
 */
public class main extends Thread {

    //各个矩阵初始化
    public static int x[][] = null;//采用分步的模式开辟数组空间
    public static int y[][] = null;//
    public static int x_m = 1;//n表示 行数 m表示 列数
    public static int x_n = 1;
    public static int z[][] = null;
    public static int M[][] = null;//计算所需的第二个矩阵M
    public static Element m[] = null;
    public static Element SK_s[] = null;
    public static Element SK_r[] = null;
    public static boolean flag = true;


    /**
     * 密钥对的初始化
     */
    public static void SKGenerator(Field Zr) {
        for (int i = 0; i < x_m; i++) {
            SK_s[i] = Zr.newRandomElement().getImmutable();
        }
        System.out.println("\n");
        for (int i = 0; i < x_n; i++) {
            SK_r[i] = Zr.newRandomElement().getImmutable();
        }
    }

    /**
     * 多线程计算VK——ProbGen
     */
    public static void calculateVK(Element[] egg_3_right_2, Element[] PK2) {
        CountDownLatch countDownLatch = new CountDownLatch(x_m);
        ExecutorService pool = Executors.newFixedThreadPool(8);
        long start = System.currentTimeMillis();
        for (int i = 0; i < x_m; i++) {
            final int d = i;
            Runnable run = () -> {
                try {
                    //模拟耗时操作
                    egg_3_right_2[d] = PK2[0].duplicate().pow(BigInteger.valueOf(x[0][d]));
                    for (int k = 1; k < x_n; k++) {
                        //为了保险起见，防止Element在运算的过程中修改了Element原本的数值，可以使用Element.duplicate()方法。
                        egg_3_right_2[d] = egg_3_right_2[d].mul(PK2[k].duplicate().pow(BigInteger.valueOf(x[k][d])));
                    }
                    countDownLatch.countDown();
                } catch (Exception e) {
                }
            };
            pool.execute(run);
        }
        try {
            countDownLatch.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        long end = System.currentTimeMillis();
        pool.shutdown();
        System.out.println("ProbGen:" + (end - start) + "ms");
    }

    /**
     * 公钥PK1 初始化
     */
    public static void calculatePK1(Element[] PK1, Element g) {
        for (int i = 0; i < x_m; i++) {
            PK1[i] = g.powZn(SK_s[i]);
        }
    }

    /**
     * 公钥PK2 初始化
     */
    public static void calculatePK2(Element[] PK2, Element g, Element h_b, Element a, Pairing bp) {
        CountDownLatch countDownLatch = new CountDownLatch(x_n);
        ExecutorService pool = Executors.newFixedThreadPool(12);
        long start = System.currentTimeMillis();
        for (int i = 0; i < x_n; i++) {
            final int d = i;
            Runnable run = () -> {
                try {
                    Element temp = g.powZn(SK_r[d]).powZn(a);
                    PK2[d] = bp.pairing(temp, h_b);
                    System.out.println("PK2:" + d + "/" + x_n);
                    countDownLatch.countDown();
                } catch (Exception e) {
                }
            };
            pool.execute(run);
        }
        try {
            countDownLatch.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        long end = System.currentTimeMillis();
        pool.shutdown();
        System.out.println("PK2:" + (end - start) + "ms");
    }

    /**
     * 计算m = s · M
     */
    public static void calculate_m() {
        for (int i = 0; i < x_n; i++) {
            Element sum = SK_s[0].mul(M[0][i]);
            for (int j = 1; j < x_m; j++) {
                sum = sum.add(SK_s[j].mul(M[j][i]));
            }
            System.out.println("m=s·M" + i + "/" + x_n);
            m[i] = sum;
        }
    }

    /**
     * 计算π
     */
    public static void calculate_pai(Element[] pai, Element delta, Element g) {
        CountDownLatch countDownLatch = new CountDownLatch(x_m);
        ExecutorService pool = Executors.newFixedThreadPool(8);
        long start = System.currentTimeMillis();
        for (int i = 0; i < x_m; i++) {
            final int j = i;
            Runnable run = () -> {
                try {
                    Element temp1 = delta.mul(m[0]).mul(x[0][j]);
                    Element temp2 = SK_r[0].mul(x[0][j]);
                    Element sum_g;
                    for (int k = 1; k < x_n; k++) {
                        temp1 = temp1.add(delta.mul(m[k]).mul(x[k][j]));
                        temp2 = temp2.add(SK_r[k].mul(x[k][j]));
                    }
                    sum_g = temp1.add(temp2);
                    pai[j] = g.powZn(sum_g);
                    countDownLatch.countDown();
                } catch (Exception e) {
                }
            };
            pool.execute(run);
        }
        try {
            countDownLatch.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        long end = System.currentTimeMillis();
        pool.shutdown();
        System.out.println("Compute:" + (end - start) + "ms");
    }

    /**
     * 计算左半边egg
     */
    public static void calculate_egg_pai_j_h(Element[] egg_pai_j_h, Element[] pai, Element a, Pairing bp, Element h_b) {
        CountDownLatch countDownLatch = new CountDownLatch(x_m);
        ExecutorService pool = Executors.newFixedThreadPool(8);
        long start = System.currentTimeMillis();
        for (int i = 0; i < x_m; i++) {
            final int j = i;
            Runnable run = () -> {
                try {
                    egg_pai_j_h[j] = pai[j].powZn(a);
                    egg_pai_j_h[j] = bp.pairing(egg_pai_j_h[j], h_b);
                    System.out.println("egg_pai_j_h:" + j + "/" + x_m);
                    countDownLatch.countDown();
                } catch (Exception e) {
                }
            };
            pool.execute(run);
        }
        try {
            countDownLatch.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        long end = System.currentTimeMillis();
        pool.shutdown();
        System.out.println("egg_pai_j_h:" + (end - start) + "ms");
    }

    /**
     * 计算右半边PK1累乘+乘方
     */
    public static void calculate_sum_PK1_y(Element[] sum_PK1_y, Element[] PK1) {
        CountDownLatch countDownLatch = new CountDownLatch(x_m);

        ExecutorService pool2 = Executors.newFixedThreadPool(4);//【最多8个，再多温度太高】
        for (int i = 0; i < x_m; i++) {
            final int j = i;
            Runnable run = () -> {
                try {
                    //开始等式右半部分
                    //先算累乘PK1^y
                    sum_PK1_y[j] = PK1[0].pow(BigInteger.valueOf(y[0][j]));//ELement必须初始化且不能自定义
                    for (int k = 1; k < x_m; k++) {
                        //TODO
                        sum_PK1_y[j] = sum_PK1_y[j].mul(PK1[k].duplicate().pow(BigInteger.valueOf(y[k][j])));//累乘
                    }
                    System.out.println("sum_PK1_y:" + j + "/" + x_m);
                    //CountDownLatch类计数减一
                    //注意，CountDownLatch类的实例要在新建多线程之前，然后入参需要
                    //统一完成后才往下执行的线程数。然后每个线程执行完后，或者
                    //部分执行完后，调用之前创建的CountDownLatch类的实例的countDown方法
                    countDownLatch.countDown();
                } catch (Exception e) {
                }
            };
            pool2.execute(run);
        }
        pool2.shutdown();

        try {
            //下面这句代码，CountDownLatch就阻塞在这里了
            //直到countDown()到0了（从构造入参的线程数开始减）
            //也即是所有线程都countDown了，
            //则解除阻塞，代码继续往下执行
            //注意，这句代码放在多线程countDown之后，多线程全部完成后
            //继续往下执行的代码之前，起一个分界线的作用
            countDownLatch.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    /**
     * 主函数入口
     */
    @SuppressWarnings("unchecked")
    public static void main(String[] args) {
        while (true) {
//            Pairing bp = PairingFactory.getPairing("a.properties");
            PairingFactory.getInstance().setUsePBCWhenPossible(true);
            int rBits = 16;//160
            int qBits = 128;//512
            TypeACurveGenerator pg = new TypeACurveGenerator(rBits, qBits);
            PairingParameters pp = pg.generate();

            Pairing bp = PairingFactory.getPairing(pp);

            Field G1 = bp.getG1();
            Field G2 = bp.getG2();
            Field Zr = bp.getZr();


            inputXandZ();
            if (isNull()) {
                System.out.println("输入矩阵不能为空");
                continue;
            }
            calculateNewX();
            calculateY();

            Element g = G1.newRandomElement().getImmutable();
            Element h = G2.newRandomElement().getImmutable();
            Element delta = Zr.newRandomElement().getImmutable();
            Element a = Zr.newRandomElement().getImmutable();
            Element b = Zr.newRandomElement().getImmutable();
            Element h_b = h.powZn(b).getImmutable();
            Element h_delta = h.powZn(delta).getImmutable();
            Element h_delta_b = h_delta.powZn(b).getImmutable();

            //密钥对的初始化
            SKGenerator(Zr);

            Element egg_3_right_2[] = new Element[x_m];
            Element PK2[] = new Element[x_n];
            Element PK1[] = new Element[x_m];
            Element pai[] = new Element[x_m];
            Element egg_pai_j_h[] = new Element[x_m];
            Element sum_PK1_y[] = new Element[x_m];


            //计算PK2
            calculatePK2(PK2, g, h_b, a, bp);

            //多线程计算VK
            calculateVK(egg_3_right_2, PK2);

            //公钥PK1 初始化
            calculatePK1(PK1, g);

            //m = s·M
            calculate_m();

            //计算π
            calculate_pai(pai, delta, g);

            //计算egg_pai_j_h双线性对
            calculate_egg_pai_j_h(egg_pai_j_h, pai, a, bp, h_b);

            //计算右半边sum_PK1_y---PK1部分
            calculate_sum_PK1_y(sum_PK1_y, PK1);

            //计算右半边最终结果并判断是否相等
            CountDownLatch countDownLatch = new CountDownLatch(x_m);
            ExecutorService pool = Executors.newFixedThreadPool(17);//TODO PK1这块多线程还需要再优化优化
            long start5 = System.currentTimeMillis();
            for (int i = 0; i < x_m; i++) {
                final int j = i;
                Runnable run = () -> {
                    try {
                        Element sum_PK1_y_a = sum_PK1_y[j].powZn(a);//m次幂运算
                        Element egg_PK1_y_h_delta = bp.pairing(sum_PK1_y_a, h_delta_b);//【别开太多线程】
                        Element egg_ab_PKparthdelta_vk = egg_PK1_y_h_delta.mul(egg_3_right_2[j]);
                        if (egg_pai_j_h[j].isEqual(egg_ab_PKparthdelta_vk)) {
                        } else {
                            flag = false;
                            System.out.println("No" + j);
                        }
                        countDownLatch.countDown();
                    } catch (Exception e) {
                    }
                };
                pool.execute(run);
            }
            pool.shutdown();

            try {
                countDownLatch.await();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            long end5 = System.currentTimeMillis();
            System.out.println("Verify:" + (end5 - start5) + "ms" + Runtime.getRuntime().availableProcessors());
            System.out.println(flag);
            flag = true;

        }
    }

    /**
     * 重写run方法
     */
    @Override
    public void run() {
        for (int i = 0; i < x_m; i++) {
            System.out.println(i);
        }
    }

    /**
     * 矩阵乘法Y尖 = M * X尖
     */
    public static void calculateY() {
        CountDownLatch countDownLatch = new CountDownLatch(x_m);
        ExecutorService pool = Executors.newFixedThreadPool(17);
        for (int i = 0; i < x_m; i++) {
            final int d = i;
            Runnable run = () -> {
                try {
                    for (int j = 0; j < x_m; j++) {
                        y[d][j] = 0;
                        for (int k = 0; k < x_n; k++) {
                            y[d][j] += M[d][k] * x[k][j];
                        }
                    }
                    System.out.println(d + "/" + x_m);
                    countDownLatch.countDown();
                } catch (Exception e) {
                }
            };
            pool.execute(run);
        }
        pool.shutdown();

        try {
            countDownLatch.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    //TODO 矩阵输入-----M
    public static void inputXandZ() {
        Scanner sc = new Scanner(System.in);
        System.out.println("输入X矩阵列的行数n：");//m表示 行数 n表示 列数
        x_n = sc.nextInt();//   ver3开始启用
        System.out.println("输入X矩阵列的列数m：");//m表示 行数 n表示 列数
        x_m = sc.nextInt();
        x = new int[x_n][x_m];
        z = new int[x_n][x_m];
        y = new int[x_m][x_m];
        //初始化SK密钥对
        SK_s = new Element[x_m];//s 对应 m
        SK_r = new Element[x_n];//r 对应 n
        //初始化m
        m = new Element[x_n];

        //TODO 初始化矩阵M
        M = new int[x_m][x_n];
        for (int i = 0; i < x_m; i++) {
            for (int j = 0; j < x_n; j++) {
                M[i][j] = (int) (1 + Math.random() * (10 - 1 + 1));
            }
        }

        for (int j = 0; j < x_n; j++) {//行
            for (int i = 0; i < x_m; i++) {
//                System.out.println("输入X矩阵第 " + (j + 1) + " 行的第 " + (i + 1) + " 个值");
//                x[j][i] = sc.nextInt();
                x[j][i] = (int) (1 + Math.random() * (10 - 1 + 1));//(数据类型)(最小值+Math.random()*(最大值-最小值+1))
            }
        }
        for (int j = 0; j < x_n; j++) {//行
            for (int i = 0; i < x_m; i++) {
//                System.out.println("输入Z矩阵第 " + (j + 1) + " 行的第 " + (i + 1) + " 个值");
//                z[j][i] = sc.nextInt();
                z[j][i] = 1;
            }
        }
        for (int i = 0; i < x_m; i++) {
            for (int j = 0; j < x_m; j++) {
                y[i][j] = 0;
            }
        }
    }

    //X尖计算
    public static void calculateNewX() {
        System.out.println("X尖 = X + Z:");
        for (int j = 0; j < x_n; j++) {
            for (int i = 0; i < x_m; i++) {
                x[j][i] = x[j][i] + z[j][i];
//                System.out.print(x[j][i] + "\t");
            }
//            System.out.print("\n");
        }
    }

    //判空
    public static boolean isNull() {
        if (x_n == 0 || x_m == 0) {
            return true;
        } else {
            return false;
        }
    }
}

