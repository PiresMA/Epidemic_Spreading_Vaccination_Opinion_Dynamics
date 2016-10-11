/** Authors: M.A.Pires and N.Crokidakis
*** Program to simulate the model presented in paper "Dynamics of epidemic spreading with vaccination: impact of social pressure and engagement"
*** Specifically this code yields the results shown in Fig.3-a, but with small modifications all the other results presented in our paper can be obtained.
*** Do not hesitate to contact us: pires.ma.fisica@gmail.com ; nuno@if.uff.br

epidemic states:
*** Vaccinated  (V): +1
*** Susceptible (S): 0
*** Infected    (I): -1

Opinion states:
*** opinion A,  pro-vaccine: +1
*** opinion B, anti-vaccine: -1

output: one file with a matrix with:
*** column 1: all values of lambda
*** column 2: S_ave(lambda) average over samples and time in the steady state
*** column 3: I_ave(lambda) average over samples and time in the steady state
*** column 4: V_ave(lambda) average over samples and time in the steady state
*** column 5: A_ave(lambda) average over samples and time in the steady state
**/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define N        10000
#define tmax       500
#define t_steady   400 // Above this instant we start to compute the average over time.
#define nsamples    10

#define so 0.99        // initial density of Susceptible agents
#define io 0.01        // initial density of Infected agents
#define vo 0.00        // initial density of Vaccinated agents

#define pao 0.60       // initial density of pro-vaccine agents. In the paper we used pao=D.

#define alpha 0.1
#define phi   0.1
#define gamma 0.2

#define dlam 0.05

#define simu 1         // label of the simulation

int main()
{
    int op[N],   op_temp[N];
    int state[N],state_temp[N];

    int    A_cont, B_cont;
    double A_cum,  B_cum;
    double A_ave,  B_ave;
    double A[tmax],B[tmax];

    int    S_cont, I_cont, V_cont;
    double S_cum,  I_cum,  V_cum;
    double S_ave,  I_ave,  V_ave;
    double S[tmax],I[tmax],V[tmax];

    // more auxiliary variables
    double lambda;

    int soma_op;
    int n1,cont,cont2,count;
    int node1,node2,node3;
    int i,j,l,t,sample;

    double aux3,aux2,aux;
    double auxA,auxB;
    double r,num;


    FILE *fp,*fq;
    char out_name[100],out_name2[100];

    //sprintf(out_name,"op_ep_time_tmax%d_N%d_pao%3.2lf_alpha%3.2lf_gamma%3.2lf_lambda%3.2lf_phi%3.2lf.dat",tmax, N,pao,alpha,gamma,lambda,phi);
    //fp=fopen(out_name,"w");

  sprintf(out_name2,"op_ep_phtr_simu%d_tmax%d_N%d_pao%3.2lf_alpha%3.2lf_gamma%3.2lf_phi%3.2lf.dat",simu,tmax, N,pao,alpha,gamma,phi);
  fq=fopen(out_name2,"w");

    for( lambda = 0 ; lambda<=1.0 ; lambda+=dlam  )
    {
        printf("lambda %3.2lf\n",lambda);

        srand(time(NULL));

        for(t=0; t<tmax; t++)
        {
            S[t]=0.0;
            I[t]=0.0;
            V[t]=0.0;
            A[t]=0.0;
            B[t]=0.0;
        }


        S_ave=0.0;
        I_ave=0.0;
        V_ave=0.0;
        A_ave=0.0;
        B_ave=0.0;

        for(sample=0; sample<nsamples; sample++)
        {

            S_cont=0;
            I_cont=0;
            V_cont=0;
            A_cont=0;
            B_cont=0;

            S_cum=0.0;
            I_cum=0.0;
            V_cum=0.0;
            A_cum=0.0;
            B_cum=0.0;

            count=0;

            aux =0.0;
            aux2=0.0;
            aux3=0.0;
            auxA=0.0;
            auxB=0.0;

            /*****    Initialization ***/
            // 1  = Vacinatted;  0  = Susceptible;  -1 = Infected

            for(i=0; i<N; i++)
            {
                r=((double) rand()/RAND_MAX);
                if(r<so)
                {
                    state[i]=0;
                    state_temp[i]=0;
                    S_cont++;
                }
                else
                {
                    state[i]=-1;
                    state_temp[i]=-1;
                    I_cont++;
                }
            } // Note that vo=0 then V_cont=0
            //printf("So=%d\tIo=%d\n",S_cont,I_cont);

            // +1  = pro-vaccine;  -1 = anti-vaccine
            for(i=0; i<N; i++)
            {
                r=((double) rand()/RAND_MAX);
                if(r<pao)
                {
                    op[i]=1;
                    op_temp[i]=1;
                    A_cont++;
                }
                else
                {
                    op[i]=-1;
                    op_temp[i]=-1;
                    B_cont++;
                }
            }
            //printf("So=%d\tIo=%d\n",S_cont,I_cont);

            S[0]+=((double) S_cont);
            I[0]+=((double) I_cont);
            V[0]+=((double) V_cont);
            A[0]+=((double) A_cont);
            B[0]+=((double) B_cont);


            /************  Dynamics ***********/

            for(t=1; t<tmax; t++)
            {
                //################################################
                /************ (I/II) Opinion Dynamics ***********/
                //################################################
                for(node1=0; node1<N; node1++)
                {
                    /************  Choose two more agents   ***********/
                    do
                    {
                        r=((double) rand()/RAND_MAX);
                        node2=((int)N*r);
                    }
                    while(node2==node1);

                    do
                    {
                        r=((double) rand()/RAND_MAX);
                        node3=((int)N*r);
                    }
                    while(node3==node1 || node3==node2 );
                    //printf("node3, node2, node1: %d, %d, %d\n",node3,node2,node1);

                    /************ Apply the majority-rule   ***********/
                    soma_op=op[node1]+op[node2]+op[node3];
                    if( abs(soma_op) < 3 )
                    {
                        if(soma_op>0)
                        {
                            op_temp[node1]=+1;
                            op_temp[node2]=+1;
                            op_temp[node3]=+1;
                        }
                        else
                        {
                            op_temp[node1]=-1;
                            op_temp[node2]=-1;
                            op_temp[node3]=-1;
                        }
                    }

                } //loop for agents

                for(j=0; j<N; j++) op[j]=op_temp[j]; // synchronous or parallel update.

                //################################################
                /************ (II/II) Epidemic Dynamics *********/
                //################################################

                for(node1=0; node1<N; node1++)
                {
                    /************ Pro-Vaccine ***********/
                    if( op[node1]==+1 )
                    {
                        if(state[node1]==0)  // if state(node1)==S
                        {
                            num=((double) rand()/RAND_MAX);
                            if(num<=gamma)
                            {
                                state_temp[node1]=+1; // S---->V with a rate gamma
                            }
                            else
                            {
                                do
                                {
                                    r=((double) rand()/RAND_MAX);
                                    node2=((int)N*r);
                                }
                                while(node2==node1);

                                if(state[node2]==-1) // if state(node2)==I
                                {
                                    num=((double) rand()/RAND_MAX);
                                    if(num<=lambda)
                                    {
                                        state_temp[node1]=-1; // S+I---->2I with a rate lambda
                                    }
                                }
                            }
                        }

                        if(state[node1]==-1) // if state(node1)==I
                        {
                            num=((double) rand()/RAND_MAX);
                            if(num<=alpha)
                            {
                                state_temp[node1]=0; // I---->S with rate alpha
                            }
                        }

                        if(state[node1]==1)  // if state(node1)==V
                        {
                            num=((double) rand()/RAND_MAX);
                            if(num<=phi)
                            {
                                state_temp[node1]=0; // V---->S with rate phi
                            }
                        }

                    }

                    /************ Anti-Vaccine ***********/
                    if( op[node1]==-1 )
                    {
                        if(state[node1]==0)  // if state(node1)==S
                        {
                            do
                            {
                                r=((double) rand()/RAND_MAX);
                                node2=((int)N*r);
                            }
                            while(node2==node1);

                            if(state[node2]==-1) // if state(node2)==S
                            {
                                num=((double) rand()/RAND_MAX);
                                if(num<=lambda)
                                {
                                    state_temp[node1]=-1; // S+I---->2I with a rate lambda
                                }
                            }
                        }

                        if(state[node1]==-1) // if state(node1)==I
                        {
                            num=((double) rand()/RAND_MAX);
                            if(num<=alpha)
                            {
                                state_temp[node1]=0; // I---->S with rate alpha
                            }
                        }

                        if(state[node1]==+1) // if state(node1)==V
                        {
                            num=((double) rand()/RAND_MAX);
                            if(num<=phi)
                            {
                                state_temp[node1]=0; // V---->S with rate phi
                            }
                        }
                    }

                } //loop for agents


                for(j=0; j<N; j++) state[j]=state_temp[j]; // synchronous or parallel update.

                S_cont=0;
                I_cont=0;
                V_cont=0;
                A_cont=0;
                B_cont=0;


                for(i=0; i<N; i++)
                {
                    if(state[i]==1 ) V_cont++;
                    if(state[i]==0 ) S_cont++;
                    if(state[i]==-1) I_cont++;
                    if(op[i]==+1) A_cont++;
                    if(op[i]==-1) B_cont++;
                }

                S[t]+=((double) S_cont);
                I[t]+=((double) I_cont);
                V[t]+=((double) V_cont);
                A[t]+=((double) A_cont);
                B[t]+=((double) B_cont);

                if(t>t_steady)
                {
                    S_cum+=S_cont;
                    I_cum+=I_cont;
                    V_cum+=V_cont;
                    A_cum+=A_cont;
                    B_cum+=B_cont;
                    count++;
                }

            } //loop for t

            // Averages over time (t>t_steady):
            aux =((double) S_cum/(N*count));
            aux2=((double) I_cum/(N*count));
            aux3=((double) V_cum/(N*count));
            auxA=((double) A_cum/(N*count));
            auxB=((double) B_cum/(N*count));

            // To compute the average over samples
            S_ave+=aux;
            I_ave+=aux2;
            V_ave+=aux3;
            A_ave+=auxA;
            B_ave+=auxB;

        } //loop for samples


//  for(t=0;t<tmax;t++)
//   {
//     S[t]/=N*nsamples;
//     I[t]/=N*nsamples;
//     V[t]/=N*nsamples;
//     A[t]/=N*nsamples;
//     B[t]/=N*nsamples;
//     fprintf(fp,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\n",t,S[t],I[t],V[t],A[t],B[t]);
//    }
//  fclose(fp);

 fprintf(fq,"%8.6lf\t%8.6lf\t%8.6lf\t%8.6lf\t%8.6lf\n",lambda,(S_ave/nsamples),(I_ave/nsamples),(V_ave/nsamples),(A_ave/nsamples));

    } // End of lambda loop

    fclose(fq);

    return(0);
}


