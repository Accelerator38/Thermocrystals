#include "particle_distribution.h"
#include "ui_particle_distribution.h"
#include <iostream>
#include <string>
#include <stdio.h>
#include <ctime>
#include <math.h>
#include <random>
#include <cassert>
#include <chrono>
#include <algorithm>
#include <thread>
#include <mutex>
#include <iostream>
#include <fstream>

using namespace std;

class RandomGenerator
{
public:
    static std::mt19937 & getMt19937();

private:
    RandomGenerator();
    ~RandomGenerator() {}
    static RandomGenerator& instance();
    RandomGenerator(RandomGenerator const&) = delete;
    RandomGenerator& operator= (RandomGenerator const&) = delete;
    std::mt19937 mMt;
};
RandomGenerator::RandomGenerator() {
    std::random_device rd;
    if (rd.entropy() != 0) {
        std::seed_seq seed{ rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd() };
        mMt.seed(seed);
    }
    else {
        auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        mMt.seed(seed);
    }
}
RandomGenerator& RandomGenerator::instance() {
    static RandomGenerator s;
    return s;
}
std::mt19937 & RandomGenerator::getMt19937() {
    return RandomGenerator::instance().mMt;
}
Particle_distribution::Particle_distribution(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::Particle_distribution)
{
    ui->setupUi(this);
}
Particle_distribution::~Particle_distribution()
{
    delete ui;
}
double Particle_distribution::random_generator()
    {
    if (ui->radio_Linear->isChecked())
        {
        return (2.0*double(rand())/RAND_MAX-1) ;
        }
    if (ui->radio_Gauss->isChecked())
        {
        double sigma = ui->sigma->text().toInt();
        double MO = ui->MO->text().toInt();
        std::normal_distribution<double> distribution(MO,sigma);
        return distribution(generator) ;
        }
    if (ui->radio_Root->isChecked())
        {
        return sin(3.14159*(double(rand())/RAND_MAX-0.5)) ;
        }
    if (ui->radio_Binary->isChecked())
        {
        int random_array[2] = {-1, 1} ;
        return random_array[rand()%2] ;
        }
    return 0 ;
    }


void Particle_distribution::on_Start_Wall_clicked()
{
    int grid = ui->grid_size->text().toInt();
    int x_size = ui->grid_size->text().toInt();
    double t_max = ui->Time_wall->text().toDouble();
    double dt = 1;
    double l = ui->lenght->text().toDouble();
    double Nw = ui->Walls_number->text().toInt();
    double I_prob = ui->Interaction_probability->text().toDouble();
    double c= ui->Speed->text().toDouble();
    double G=ui->Gravity->text().toDouble();
    int N = ui->particle_number->text().toInt();

    int time_steps=round(t_max/dt) ;

    QVector<double> x_1(2*grid), y_1(2*grid), x_2(2*grid), y_2(2*grid), x_3(t_max/dt), y_3(t_max/dt) ;
    std::vector<double> S(N), V(N), W(2*grid) ;
    std::vector<double> Ssumm (N), Vsumm (N), R (t_max/dt), y_1_plus (N)  ;
    std::vector<std::vector<double>> P(N, std::vector<double>(time_steps)); //матрица Плотности вещества

    int k=0 ;


    ui->widget->addGraph();
    ui->widget->graph(0)->setData(x_1, y_1);
    ui->widget->graph(0)->setPen(QColor(50, 50, 50, 255));
    ui->widget->graph(0)->setLineStyle(QCPGraph::lsNone);
    ui->widget->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    ui->widget->xAxis->setLabel("X, координата");
    ui->widget->yAxis->setLabel("N, количество частиц");
    ui->widget->xAxis->setRange(-x_size, x_size);
    ui->widget->yAxis->setRange(0, 0.05 );
    ui->widget->replot();

    ui->widget_2->addGraph();
    ui->widget_2->graph(0)->setData(x_2, y_2);
    ui->widget_2->xAxis->setLabel("X, координата");
    ui->widget_2->yAxis->setLabel("N, отн. ед.");
    ui->widget_2->xAxis->setRange(-x_size, x_size);
    ui->widget_2->yAxis->setRange(-1.1,1.1 );
    ui->widget_2->replot();

    std::mt19937 &mt = RandomGenerator::getMt19937();
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::mt19937 &mt2 = RandomGenerator::getMt19937();
    std::uniform_real_distribution<double> dist2(0.0, 1.0);

for (int i=0; i<2*grid; i++) //барьеры
    {
    if (i<=Nw)
        {
        W[i]=-1 ;
        }
    else
        {
        W[i]=1 ;
        }
    }


    for (int i=0 ; i<N; i++) // начальные условия
        {
        V[i]=c*random_generator() ;
        S[i]=l*(2.0*double(rand())/RAND_MAX-1) ;
        }




for (double t=0; t<=t_max; t+=dt) //основной цикл по времени
    {
    srand(unsigned(time(nullptr)));
    random_shuffle(W.begin(), W.end());

    for (int i=0 ; i<N; i++) //координаты
        {
        S[i]=S[i]+V[i]*dt/c ;
        }

    for (int i=0 ; i<N; i++) //скорости
        {
        if (dist2(mt2)*100<I_prob)
            {
            V[i]=V[i]*W[int(round(S[i])) + grid] ;
            }
        }
/*
    for (int i=0 ; i<N; i++) //Притяжение
        {
        if (S[i]<0)
            {
            V[i]=V[i]+G ;
            }
        else
            {
            V[i]=V[i]-G ;
            }
        }
*/
        for (int i=0 ; i<2*grid; i++) //визуализация
            {
            x_1[i]=i-grid;
            y_1[i]=0;
            }

        for (int i=0 ; i<N; i++)
            {
            int X_grid=int(round(S[i]))+grid;
            if (X_grid>=0 && X_grid<2*grid)
                {
                y_1[X_grid]++;
                }
            }

        for (int i=0 ; i<2*grid; i++) //визуализация
            {
            y_1[i]=y_1[i]/N;
            P[i][k]=y_1[i];
            }

        for (int i=0 ; i<2*grid; i++) //представление результатов
            {
            R[k]+=S[i]*S[i];
            }

        k++ ;

        ui->widget->graph(0)->setData(x_1, y_1);
        ui->widget->replot();
    }



    for (int i=0; i<2*grid; i++)
        {
        y_1[i] = y_1[i]/N ;
        }

    for (int i=0 ; i<t_max/dt; i+=dt) //представление результатов
        {
        y_3[i]=R[i] ;
        x_3[i]=i ;
        }

    double y_3_max=0 ;
    for (int i=0; i<t_max/dt; i++) //поиск максимума
        {
        if (y_3_max<y_3[i])
            {
                y_3_max=y_3[i] ;
            }
        }
    for (int i=0; i<t_max/dt; i++)
        {
        y_3[i]=y_3[i]/y_3_max ;
        }

ui->widget_2->addGraph();
ui->widget_2->graph(0)->setData(x_1, y_1);
ui->widget_2->xAxis->setLabel("X, координата");
ui->widget_2->yAxis->setLabel("N, отн. ед.");
ui->widget_2->xAxis->setRange(-x_size, x_size);
ui->widget_2->yAxis->setRange(0,1.1 );
ui->widget_2->replot();




std::ofstream out;
double tt=0;
out.open("D:\\PMatrix.txt"); // окрываем файл для записи
if (out.is_open())
    {
    for (int i=0; i<time_steps; i++)
        {
        if (i>=tt)
            {
            tt+=time_steps/10 ;
            for (int j=0; j<N; j++)
                {
                out << P[j][i] <<'\t' ;
                }
            out << std::endl ;
            }
        }
    }
out.close();
}

void Particle_distribution::on_Start_Chains_clicked()
{
    int N=ui->particle_number_chain->text().toInt();
    int N_chains=ui->Number_chains->text().toInt();
    double k=ui->Integration_Step->text().toDouble();
    double E=ui->Epsilon->text().toDouble();
    double Signal=ui->Signal_Level->text().toDouble();
    double time1=ui->Time1->text().toDouble(); //время стабилизации температур
    double time2=ui->Time2->text().toDouble(); //время включения разрывов
    double t_max=ui->Time_Chain->text().toDouble(); // время процесса
    double L=ui->lenght_2->text().toDouble();
    double y=ui->y_speed->text().toDouble();
    double C=100 ;
    double a=1 ;
    double m=1 ;
    double L_noise=50 ;
    double we=sqrt(C/m) ;
    double b=a*(1+E) ;
    double c=a*we;
    double Te=2*3.14159/we ;
    double dt=k*Te ;
    double V_signal=1+Signal ;
    double V_noise=1;
    double Vrand=y*c;
    int T_number=round(t_max/dt) ;

    QVector<double> x_1(N), x_2(T_number), y_1(N), y_2(T_number), y_3(T_number),y_4 (T_number) ;
    std::vector<double> S_summ (T_number), V_summ (T_number), T_summ(N), T_local(N), R_summ(T_number), W_summ(N) ;
    std::vector<std::vector<double>> T(N, std::vector<double>(T_number)); //матрица Температур
    std::mt19937 &mt = RandomGenerator::getMt19937();
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    int Threads_Number=ui->Threads_number->text().toDouble(); //число потоков
    std::mutex BLOСK ; //блокирование потоков

    auto calc_chains = [&]() -> void
        {
        std::vector<double> S(N), V(N),  W(N), S_0s(2*L+1), V_0s(2*L+1) ;
        for (int j=0; j<N_chains/Threads_Number; j++)
            {
            int VV_S0=0 ;
            std::vector<std::vector<double>> R(N, std::vector<double>(T_number)); //матрица радиусов

            for (int k=0 ; k<N; k++) //задаем начальные условия (t=0)
                {
                S[k]=0 ;
                V[k]=V_noise*Vrand*(dist(mt)-0.5) ;
                }                       
            for (int k=0 ; k<=2*L; k++)
                {
                S_0s[k]=0 ;
                V_0s[k]=V_signal*Vrand*(dist(mt)-0.5) ;
                }

            int time_count=0 ;
            for (double t=dt; t<=t_max; t+=dt) //цикл по времени
                {

                if (t<time1) // начало выравнивания по кинетической и потенциальной энергии
                    {
                    for (int k=1 ; k<N-1; k++)
                        {
                        V[k]+=dt*we*we*((S[k-1]-S[k])+(S[k+1]-S[k])) ;
                        }
                    V[0]+=dt*we*we*((S[1]-S[0])+(S[N-2]-S[0])) ;
                    V[N-1]=V[0] ; //граничные условия
                    for (int k=0 ; k<N; k++)
                        {
                        S[k]+=dt*V[k] ;
                        }

                    for (int k=1 ; k<2*L-1; k++)
                        {
                        V_0s[k]+=dt*we*we*((S_0s[k-1]-S_0s[k])+(S_0s[k+1]-S_0s[k])) ;
                        }
                    V_0s[0]+=dt*we*we*((S_0s[1]-S_0s[0])+(S_0s[2*L-2]-S[0])) ;
                    V_0s[2*L-1]=V_0s[0] ; //граничные условия
                    for (int k=0 ; k<2*L; k++)
                        {
                        S_0s[k]+=dt*V_0s[k] ;
                        }
                    } //конец выравнивания энергии


                if ((t>=time1)&(VV_S0==0)) //задание сигнала
                    {
                    VV_S0=1;
                    /*
                    for (int k=0 ; k<=2*L; k++)
                        {
                        S[N/2-L+k]=0 ;
                        V[N/2-L+k]=V_signal*Vrand*(dist(mt)-0.5) ;
                        }
                    */
                    for (int k=1 ; k<=2*L; k++)
                        {
                        S[N/2-L+k]=S_0s[k] ;
                        V[N/2-L+k]=V_0s[k] ;
                        }

                    }


                if ((t>=time1)&(t<time2)) // начало теплопереноса без разрывов
                    {
                    for (int k=1 ; k<N-1; k++)
                        {
                        V[k]+=dt*we*we*((S[k-1]-S[k])+(S[k+1]-S[k])) ;
                        }
                    V[0]+=dt*we*we*((S[1]-S[0])+(S[N-2]-S[0])) ;
                    V[N-1]=V[0] ; //граничные условия

                    for (int k=0 ; k<N; k++)
                        {
                        S[k]+=dt*V[k] ;
                        }

                    double Noise=0;
                    for (int k=0 ; k<N; k++) //находим уровень теплового шума
                        {
                        if ((k<L_noise)||(k>N-L_noise))
                            {
                            Noise+=V[k]*V[k];
                            }
                        }
                    Noise=Noise/(L_noise*2) ;

                    for (int k=0 ; k<N; k++)
                        {
                        R[k][time_count]=(V[k]*V[k]-Noise)*(N/2-k)*(N/2-k) ;
                        }
                    for (int k=0 ; k<N; k++)
                        {
                        V_summ[time_count]+=V[k]*V[k];
                        }
                    } //конец второго периода


                if (t>=time2) // начало теплопереноса с разрывами
                    {
                    for (int k=0 ; k<N-2; k++)
                        {
                        if (abs(S[k+1]-S[k])>E)
                            {
                            W[k]=0 ;
                            W_summ[k]++ ;
                            }
                        else
                            {
                            W[k]=1 ;
                            }
                        }

                    if (abs(S[N-2]-S[0])>E)
                        {
                        W[N-2]=0 ;
                        W_summ[k]++ ;
                        }
                    else
                        {
                        W[N-2]=1 ;
                        }

                    for (int k=1 ; k<N-1; k++)
                        {
                        V[k]+=dt*we*we*(W[k-1]*(S[k-1]-S[k])+W[k]*(S[k+1]-S[k])) ;
                        }
                    V[0]+=dt*we*we*(W[0]*(S[1]-S[0])+W[N-2]*(S[N-2]-S[0])) ;
                    V[N-1]=V[0] ; //граничные условия

                    for (int k=0 ; k<N; k++)
                        {
                        S[k]+=dt*V[k] ;
                        }

                    double Noise=0;
                    for (int k=0 ; k<N; k++) //находим уровень теплового шума
                        {
                        if ((k<L_noise)||(k>N-L_noise))
                            {
                            Noise+=V[k]*V[k];
                            }
                        }
                    Noise=Noise/(2*L_noise);

                    for (int k=0 ; k<N; k++)
                        {
                        R[k][time_count]=(V[k]*V[k]-Noise)*(N/2-k)*(N/2-k) ;
                        }
                    for (int k=0 ; k<N; k++)
                        {
                        T[k][time_count]+=V[k]*V[k] ;
                        }

                    for (int k=0 ; k<N; k++)
                        {
                        V_summ[time_count]+=V[k]*V[k];
                        }
                    } //конец третьего периода


                time_count++ ;
                } //конец цикла по времени


            std::lock_guard<std::mutex> guard(BLOСK); // блокирование потоков и разблокирование при выходи из БЛОКА
            //BLOСK.lock();
            for (int k=0 ; k<N; k++)
                {

                double Noise=0;
                for (int k=0 ; k<L_noise; k++) //находим уровень теплового шума
                    {
                    Noise+=V[k]*V[k];
                    }
                for (int k=N-L_noise ; k<N ; k++)
                    {
                    Noise+=V[k]*V[k];
                    }
                Noise=Noise/(L_noise*2) ;


                T_summ[k]+=V[k]*V[k]-Noise ;
                }


            for (int j=0 ; j<T_number; j++)
                {
                for (int i=0 ; i<N; i++)
                    {
                    R_summ[j]+=R[i][j] ;
                    }
                }
            //BLOСK.unlock();
            }
        };




std::vector<std::thread> Threads(Threads_Number); //создание потоков
for (int i=0 ; i<Threads.size(); i++)
    {
    Threads[i] = std::thread(calc_chains);
    }
for (int i=0 ; i<Threads.size(); i++)
    {
    Threads[i].join(); //join - соединение потоков, то есть ожидание завершения потока
    }



double y1_max=0 ;
double y4_max=0 ;
for (int i=0 ; i<N; i++)
    {
    x_1[i]=i-N/2 ;
    y_1[i]=T_summ[i]/N_chains ;
    y_4[i]=W_summ[i]/(N_chains*T_number);
    if (y_1[i]>y1_max)
        {
        y1_max=y_1[i];
        }
    if (y_4[i]>y4_max)
        {
        y4_max=y_4[i];
        }
    }

for (int i=0 ; i<N; i++)
    {
    y_1[i]=y_1[i]/y1_max;
    y_4[i]=y_4[i]/y4_max;
    }

for (int j=0 ; j<T_number; j++)
    {
    y_2[j]=R_summ[j]/N_chains ;
    y_3[j]=V_summ[j]/N_chains ;
    }
double y2_max=0 ;
double y3_max=0 ;
for (int j=0 ; j<T_number; j++)
    {
    if (y_2[j]>y2_max)
        {
        y2_max=y_2[j] ;
        }
    if (y_3[j]>y3_max)
        {
        y3_max=y_3[j] ;
        }
    }
for (int j=0 ; j<T_number; j++)
    {
    x_2[j]=dt*j ;
    y_2[j]=y_2[j]/y2_max ;
    y_3[j]=y_3[j]/y3_max ;
    }

    ui->widget_2->addGraph();
    ui->widget_2->addGraph();
    ui->widget_2->graph(1)->setLineStyle(QCPGraph::lsNone);
    ui->widget_2->graph(1)->setData(x_1, y_1);
    ui->widget_2->graph(1)->setPen(QColor(200, 0, 0, 255));
    ui->widget_2->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    ui->widget_2->xAxis->setLabel("X, Номер узла");
    ui->widget_2->yAxis->setLabel("T, температура");
    ui->widget_2->xAxis->setRange(-0.5*N, 0.5*N);
    ui->widget_2->yAxis->setRange(0.0, 1.001 );
    ui->widget_2->replot();

    ui->widget_3->addGraph();
    ui->widget_3->graph(0)->setData(x_2, y_2);
    ui->widget_3->graph(0)->setPen(QColor(50, 50, 50, 255));
    ui->widget_3->graph(0)->setLineStyle(QCPGraph::lsNone);
    ui->widget_3->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    ui->widget_3->xAxis->setLabel("t, время");
    ui->widget_3->yAxis->setLabel("R, радиус инерции");
    ui->widget_3->xAxis->setRange(0, t_max+dt);
    ui->widget_3->yAxis->setRange(0,1.1 );
    ui->widget_3->replot();

    ui->widget_4->addGraph();
    ui->widget_4->graph(0)->setData(x_2, y_3);
    ui->widget_4->graph(0)->setPen(QColor(50, 50, 50, 255));
    ui->widget_4->graph(0)->setLineStyle(QCPGraph::lsNone);
    ui->widget_4->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    ui->widget_4->xAxis->setLabel("t, время");
    ui->widget_4->yAxis->setLabel("VV");
    ui->widget_4->xAxis->setRange(0, t_max+dt);
    ui->widget_4->yAxis->setRange(0,1.1 );
    ui->widget_4->replot();

    ui->widget_5->addGraph();
    ui->widget_5->graph(0)->setData(x_1, y_4);
    ui->widget_5->graph(0)->setPen(QColor(50, 50, 50, 255));
    ui->widget_5->graph(0)->setLineStyle(QCPGraph::lsNone);
    ui->widget_5->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    ui->widget_5->xAxis->setLabel("X, Номер узла");
    ui->widget_5->yAxis->setLabel("Вероятность разрыва");
    ui->widget_5->xAxis->setRange(-0.5*N, 0.5*N);
    ui->widget_5->yAxis->setRange(0,1 );
    ui->widget_5->replot();


    std::ofstream out;
    out.open("D:\\TempLast.txt");
    if (out.is_open())
        {
        for (int i=0; i<N; i++)
            {
            out << T_summ[i]<<'\t' ;
            }
        }
    out.close();

    out.open("D:\\Walls.txt"); // окрываем файл для записи
    if (out.is_open())
        {
        for (int i=0; i<N; i++)
            {
            out << W_summ[i]<<'\t' ;
            }
        }
    out.close();

    out.open("D:\\Radius.txt"); // окрываем файл для записи
    if (out.is_open())
        {
        for (int i=0; i<T_number; i++)
            {
            out << R_summ[i]<<'\t' ;
            }
        }
    out.close();



    double tt=0;


    for (int i=0; i<T_number; i++)
        {
        if (i>=tt)
            {
            tt+=T_number/10 ;
            double SummT=0 ;

            for (int j=0; j<N; j++)
                {
                SummT+=T[j][i]  ;
                }

            for (int j=0; j<N; j++)
                {
                T[j][i]=T[j][i]/SummT ;
                }
            }
        }




    tt=0;
    out.open("D:\\TempMatrix.txt"); // окрываем файл для записи
    if (out.is_open())
        {
        for (int i=0; i<T_number; i++)
            {
            if (i>=tt)
                {
                tt+=T_number/10 ;
                for (int j=0; j<N; j++)
                    {
                    out << T[j][i] <<'\t' ;
                    }
                out << std::endl ;
                }
            }
        }
    out.close();
}
