//
//  CollinearEQ.cpp
//  
//
//  Created by MJ_Noh on 8/6/21.
//

#include "CollinearEQ.hpp"

bool CollinearCalibration(vector<D2DPOINT> &IPs, vector<D3DPOINT> &GCPs, EO &eo, CAMERA_INFO &camera, bool frame)
{
    bool check_solution = false;
    int numofpts = IPs.size();
    
    if(frame)
    {
        EO p_eo = eo;
        CAMERA_INFO p_camera = camera;
        /*
        GetInitialPCfromDLT(IPs, GCPs, eo, camera);
        printf("initial eo %f\t%f\t%f\t%f\t%f\t%f\n",eo.m_Wl*RadToDeg,eo.m_Pl*RadToDeg,eo.m_Kl,eo.m_Xl,eo.m_Yl,eo.m_Zl);
        printf("initial ca %f\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\n",
               camera.m_focalLength,camera.m_ppx,camera.m_ppy,
               camera.k1,camera.k2,camera.k3,
               camera.p1,camera.p2,
               camera.a1,camera.a2);
        */
        bool check_conver = CalibrationBundle1(IPs, GCPs, p_eo, p_camera);
        if(!check_conver)
        {
            EO pp_eo = eo;
            CAMERA_INFO pp_camera;
            pp_camera.m_focalLength = camera.m_focalLength;
            pp_camera.m_ppx = 0;
            pp_camera.m_ppy = 0;
            pp_camera.k1 = 0;
            pp_camera.k2 = 0;
            pp_camera.k3 = 0;
            pp_camera.p1 = 0;
            pp_camera.p2 = 0;
            pp_camera.a1 = 0;
            pp_camera.a2 = 0;
            
            bool check_conver2 = CalibrationBundlewithoutFL(IPs, GCPs, pp_eo, pp_camera);
            
            if(!check_conver2)
            {
                eo.m_Wl = 0;
                eo.m_Pl = 0;
                eo.m_Kl = 0;
                
                eo.m_Xl = 0;
                eo.m_Yl = 0;
                eo.m_Zl = 0;
                
                printf("No solution\n");
                check_solution = false;
            }
            else
            {
                eo = pp_eo;
                camera = pp_camera;
                check_solution = true;
            }
            
        }
        else
        {
            eo = p_eo;
            camera = p_camera;
            check_solution = true;
        }
        printf("Adjusted eo %f\t%f\t%f\t%f\t%f\t%f\n",eo.m_Wl*RadToDeg,eo.m_Pl*RadToDeg,eo.m_Kl*RadToDeg,eo.m_Xl,eo.m_Yl,eo.m_Zl);
        printf("initial ca %f\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\n",
               camera.m_focalLength,camera.m_ppx,camera.m_ppy,
               camera.k1,camera.k2,camera.k3,
               camera.p1,camera.p2,
               camera.a1,camera.a2);
    }
    else
    {
        GetInitialPCfromDLT(IPs, GCPs, eo, camera);
        printf("initial eo %f\t%f\t%f\t%f\t%f\t%f\n",eo.m_Wl*RadToDeg,eo.m_Pl*RadToDeg,eo.m_Kl,eo.m_Xl,eo.m_Yl,eo.m_Zl);
        printf("initial ca %f\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\n",
               camera.m_focalLength,camera.m_ppx,camera.m_ppy,
               camera.k1,camera.k2,camera.k3,
               camera.p1,camera.p2,
               camera.a1,camera.a2);
        /*
        CalibrationBundle1(IPs, GCPs, eo, camera);
        printf("Adjusted eo %f\t%f\t%f\t%f\t%f\t%f\n",eo.m_Wl*RadToDeg,eo.m_Pl*RadToDeg,eo.m_Kl*RadToDeg,eo.m_Xl,eo.m_Yl,eo.m_Zl);
        printf("initial ca %f\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\n",
               camera.m_focalLength,camera.m_ppx,camera.m_ppy,
               camera.k1,camera.k2,camera.k3,
               camera.p1,camera.p2,
               camera.a1,camera.a2);
        
        */
        EOEstimatefromInitial(IPs, GCPs, eo, camera);
        printf("Adjusted eo %f\t%f\t%f\t%f\t%f\t%f\n",eo.m_Wl*RadToDeg,eo.m_Pl*RadToDeg,eo.m_Kl*RadToDeg,eo.m_Xl,eo.m_Yl,eo.m_Zl);
        printf("initial ca %f\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\n",
               camera.m_focalLength,camera.m_ppx,camera.m_ppy,
               camera.k1,camera.k2,camera.k3,
               camera.p1,camera.p2,
               camera.a1,camera.a2);
        
        PPAEstimatefromEO(IPs, GCPs, eo, camera);
        printf("Adjusted eo %f\t%f\t%f\t%f\t%f\t%f\n",eo.m_Wl*RadToDeg,eo.m_Pl*RadToDeg,eo.m_Kl*RadToDeg,eo.m_Xl,eo.m_Yl,eo.m_Zl);
        printf("initial ca %f\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\n",
               camera.m_focalLength,camera.m_ppx,camera.m_ppy,
               camera.k1,camera.k2,camera.k3,
               camera.p1,camera.p2,
               camera.a1,camera.a2);
        
        RadialParamEstimate(IPs, GCPs, eo, camera);
        printf("Adjusted eo %f\t%f\t%f\t%f\t%f\t%f\n",eo.m_Wl*RadToDeg,eo.m_Pl*RadToDeg,eo.m_Kl*RadToDeg,eo.m_Xl,eo.m_Yl,eo.m_Zl);
        printf("initial ca %f\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\n",
               camera.m_focalLength,camera.m_ppx,camera.m_ppy,
               camera.k1,camera.k2,camera.k3,
               camera.p1,camera.p2,
               camera.a1,camera.a2);
        
        TangentialParamEstimate(IPs, GCPs, eo, camera);
        printf("Adjusted eo %f\t%f\t%f\t%f\t%f\t%f\n",eo.m_Wl*RadToDeg,eo.m_Pl*RadToDeg,eo.m_Kl*RadToDeg,eo.m_Xl,eo.m_Yl,eo.m_Zl);
        printf("initial ca %f\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\n",
               camera.m_focalLength,camera.m_ppx,camera.m_ppy,
               camera.k1,camera.k2,camera.k3,
               camera.p1,camera.p2,
               camera.a1,camera.a2);
        
        AffinityParamEstimate(IPs, GCPs, eo, camera);
        printf("Adjusted eo %f\t%f\t%f\t%f\t%f\t%f\n",eo.m_Wl*RadToDeg,eo.m_Pl*RadToDeg,eo.m_Kl*RadToDeg,eo.m_Xl,eo.m_Yl,eo.m_Zl);
        printf("initial ca %f\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\n",
               camera.m_focalLength,camera.m_ppx,camera.m_ppy,
               camera.k1,camera.k2,camera.k3,
               camera.p1,camera.p2,
               camera.a1,camera.a2);
        
        check_solution = true;
    }
    //exit(1);
    // */
    eo.m_Wl *= RadToDeg;
    eo.m_Pl *= RadToDeg;
    eo.m_Kl *= RadToDeg;
    
    return check_solution;
}

void GetInitialPCfromDLT(vector<D2DPOINT> &IPs, vector<D3DPOINT> &GCPs, EO &eo, CAMERA_INFO &camera)
{
    int numofpts = IPs.size();
    GMA_double *A_matrix = GMA_double_create(numofpts*2,11);
    GMA_double *L_matrix = GMA_double_create(numofpts*2,1);
    
    
    for(int i=0;i<numofpts;i++)
    {
        //printf("ID %d\tIP %f\t%f\tVCP %f\t%f\t%f\n",i,IPs[i].m_X,IPs[i].m_Y,GCPs[i].m_X,GCPs[i].m_Y,GCPs[i].m_Z);
        
        A_matrix->val[i*2][0] = GCPs[i].m_X;
        A_matrix->val[i*2][1] = GCPs[i].m_Y;
        A_matrix->val[i*2][2] = GCPs[i].m_Z;
        A_matrix->val[i*2][3] = 1;
        A_matrix->val[i*2][4] = 0;
        A_matrix->val[i*2][5] = 0;
        A_matrix->val[i*2][6] = 0;
        A_matrix->val[i*2][7] = 0;
        A_matrix->val[i*2][8]  = -(GCPs[i].m_X*IPs[i].m_X);
        A_matrix->val[i*2][9]  = -(GCPs[i].m_Y*IPs[i].m_X);
        A_matrix->val[i*2][10] = -(GCPs[i].m_Z*IPs[i].m_X);
        
        L_matrix->val[i*2][0]  = IPs[i].m_X;
        
        A_matrix->val[i*2 + 1][0] = 0;
        A_matrix->val[i*2 + 1][1] = 0;
        A_matrix->val[i*2 + 1][2] = 0;
        A_matrix->val[i*2 + 1][3] = 0;
        A_matrix->val[i*2 + 1][4] = GCPs[i].m_X;
        A_matrix->val[i*2 + 1][5] = GCPs[i].m_Y;
        A_matrix->val[i*2 + 1][6] = GCPs[i].m_Z;
        A_matrix->val[i*2 + 1][7] = 1;
        A_matrix->val[i*2 + 1][8]  = -(GCPs[i].m_X*IPs[i].m_Y);
        A_matrix->val[i*2 + 1][9]  = -(GCPs[i].m_Y*IPs[i].m_Y);
        A_matrix->val[i*2 + 1][10] = -(GCPs[i].m_Z*IPs[i].m_Y);
        
        L_matrix->val[i*2 + 1][0]  = IPs[i].m_Y;
    }
    
    //GMA_double_printf(A_matrix);
    //GMA_double_printf(L_matrix);
    GMA_double *AT_matrix   = GMA_double_create(11,numofpts*2);
    GMA_double *ATA_matrix  = GMA_double_create(11,11);
    GMA_double *ATAI_matrix = GMA_double_create(11,11);
    GMA_double *ATL_matrix  = GMA_double_create(11,1);
    GMA_double *X_matrix    = GMA_double_create(11,1);
    
    GMA_double_Tran(A_matrix,AT_matrix);
    GMA_double_mul(AT_matrix,A_matrix,ATA_matrix);
    GMA_double_inv(ATA_matrix,ATAI_matrix);
    GMA_double_mul(AT_matrix,L_matrix,ATL_matrix);
    GMA_double_mul(ATAI_matrix,ATL_matrix,X_matrix);
    
    GMA_double *AX_matrix    = GMA_double_create(numofpts*2,1);
    GMA_double *V_matrix     = GMA_double_create(numofpts*2,1);
    
    GMA_double_mul(A_matrix,X_matrix,AX_matrix);
    GMA_double_sub(AX_matrix,L_matrix,V_matrix);
    
    
    
    double var_sum = 0;
    for(int i=0;i<numofpts*2;i++)
        var_sum += (V_matrix->val[i][0])*(V_matrix->val[i][0]);
    
    double sigma = sqrt(var_sum/(numofpts*2 - 11));
    
    DLTPARAM param;
    param.A1 = X_matrix->val[0][0];
    param.B1 = X_matrix->val[1][0];
    param.C1 = X_matrix->val[2][0];
    param.D1 = X_matrix->val[3][0];
    param.A2 = X_matrix->val[4][0];
    param.B2 = X_matrix->val[5][0];
    param.C2 = X_matrix->val[6][0];
    param.D2 = X_matrix->val[7][0];
    param.A3 = X_matrix->val[8][0];
    param.B3 = X_matrix->val[9][0];
    param.C3 = X_matrix->val[10][0];
    
    //printf("DLT sigma %f\n", sigma);
    //GMA_double_printf(X_matrix);
    //exit(1);
    GMA_double_destroy(A_matrix);
    GMA_double_destroy(L_matrix);
    GMA_double_destroy(AT_matrix);
    GMA_double_destroy(ATA_matrix);
    GMA_double_destroy(ATAI_matrix);
    GMA_double_destroy(ATL_matrix);
    GMA_double_destroy(X_matrix);
    GMA_double_destroy(AX_matrix);
    GMA_double_destroy(V_matrix);
    
    double L = -sqrt(param.A3*param.A3 + param.B3*param.B3 + param.C3*param.C3);
    if(L == 0.0)
        L = 0.000000000001;
    
    double xp, yp, fx, fy;
    double r11, r12, r13;
    double r21, r22, r23;
    double r31, r32, r33;
    
    xp = (param.A1*param.A3 + param.B1*param.B3 + param.C1*param.C3) / (L*L);
    yp = (param.A2*param.A3 + param.B2*param.B3 + param.C2*param.C3) / (L*L);
        
    fx = sqrt(((param.A1*param.A1 + param.B1*param.B1 + param.C1*param.C1)/ (L*L)) - xp*xp);
    fy = sqrt(((param.A2*param.A2 + param.B2*param.B2 + param.C2*param.C2)/ (L*L)) - yp*yp);
    
    //printf("L %f\txp %f\typ %f\tfx %f\tfy %f\n",L,xp,yp,fx,fy);
    
    r31 = param.A3/L;
    r32 = param.B3/L;
    r33 = param.C3/L;

    r11 = (xp*r31 - (param.A1/L)) / fx;
    r12 = (xp*r32 - (param.B1/L)) / fx;
    r13 = (xp*r33 - (param.C1/L)) / fx;

    r21 = (yp*r31 - (param.A2/L)) / fy;
    r22 = (yp*r32 - (param.B2/L)) / fy;
    r23 = (yp*r33 - (param.C2/L)) / fy;
    
    A_matrix = GMA_double_create(3,3);
    L_matrix = GMA_double_create(3,1);
    X_matrix = GMA_double_create(3,1);
    
    A_matrix->val[0][0] = param.A1;
    A_matrix->val[0][1] = param.B1;
    A_matrix->val[0][2] = param.C1;
    
    A_matrix->val[1][0] = param.A2;
    A_matrix->val[1][1] = param.B2;
    A_matrix->val[1][2] = param.C2;
    
    A_matrix->val[2][0] = param.A3;
    A_matrix->val[2][1] = param.B3;
    A_matrix->val[2][2] = param.C3;
    
    L_matrix->val[0][0] = param.D1;
    L_matrix->val[1][0] = param.D2;
    L_matrix->val[2][0] = 1.0;
    
    GMA_double *AI_matrix = GMA_double_create(3,3);
    GMA_double *AIL_matrix = GMA_double_create(3,1);
    
    GMA_double_inv(A_matrix,AI_matrix);
    GMA_double_mul(AI_matrix,L_matrix,AIL_matrix);
    
    eo.m_Xl = - AIL_matrix->val[0][0];
    eo.m_Yl = - AIL_matrix->val[1][0];
    eo.m_Zl = - AIL_matrix->val[2][0];
    
    /*
    eo.m_Wl = atan2(-r32,r33);
    eo.m_Pl = asin(r31);
    eo.m_Kl = atan2(-r21,r11);
    
    camera.m_focalLength = (fx + fy)/2.0;
    camera.m_ppx = xp;
    camera.m_ppy = yp;
    */
    /*
    eo.m_Xl = GCPs[0].m_X;
    eo.m_Yl = GCPs[0].m_Y;
    eo.m_Zl = 490000;
    */
    GMA_double_destroy(A_matrix);
    GMA_double_destroy(L_matrix);
    GMA_double_destroy(X_matrix);
    GMA_double_destroy(AI_matrix);
    GMA_double_destroy(AIL_matrix);
    
    printf("initial eo %f\t%f\t%f\t%f\t%f\t%f\n",eo.m_Wl*RadToDeg,eo.m_Pl*RadToDeg,eo.m_Kl,eo.m_Xl,eo.m_Yl,eo.m_Zl);
    printf("initial ca %f\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\n",
           camera.m_focalLength,camera.m_ppx,camera.m_ppy,
           camera.k1,camera.k2,camera.k3,
           camera.p1,camera.p2,
           camera.a1,camera.a2);
}

void EOEstimatefromInitial(vector<D2DPOINT> &IPs, vector<D3DPOINT> &GCPs, EO &eo, CAMERA_INFO &camera)
{
    int numofpts = IPs.size();
    int unknowns = 6;
    
    GMA_double *A_matrix = GMA_double_create(2*numofpts,unknowns);
    GMA_double *L_matrix = GMA_double_create(2*numofpts,1);
    GMA_double *AT_matrix = GMA_double_create(unknowns,2*numofpts);
    GMA_double *ATA_matrix = GMA_double_create(unknowns,unknowns);
    GMA_double *ATAI_matrix = GMA_double_create(unknowns,unknowns);
    GMA_double *ATL_matrix  = GMA_double_create(unknowns,1);
    GMA_double *X_matrix    = GMA_double_create(unknowns,1);
    GMA_double *AX_matrix    = GMA_double_create(numofpts*2,1);
    GMA_double *V_matrix     = GMA_double_create(numofpts*2,1);
    
    bool check_stop = false;
    int max_iteration = 100;
    double max_correct = 0.000001;
    int iter = 0;
    while(!check_stop && iter < max_iteration)
    {
        iter++;
        
        UPARAMS uparam;
        PDCS pdcs;
        
        for(int i=0;i<numofpts;i++)
        {
            uparam.f = camera.m_focalLength;
            uparam.k1 = 0;
            uparam.k2 = 0;
            uparam.k3 = 0;
            uparam.p1 = 0;
            uparam.p2 = 0;
            uparam.xp = camera.m_ppx;
            uparam.yp = camera.m_ppy;
            uparam.a1 = 0;
            uparam.a2 = 0;
            
            uparam.X  = eo.m_Xl;
            uparam.Y  = eo.m_Yl;
            uparam.Z  = eo.m_Zl;
            uparam.omega = eo.m_Wl;
            uparam.phi   = eo.m_Pl;
            uparam.kappa = eo.m_Kl;
            
            uparam.xa = IPs[i].m_X;
            uparam.ya = IPs[i].m_Y;
            uparam.XA = GCPs[i].m_X;
            uparam.YA = GCPs[i].m_Y;
            uparam.ZA = GCPs[i].m_Z;
            
            GetPDCs(uparam,pdcs);
            
            for(int j=0;j<unknowns;j++)
            {
                double sign;
                if(j>2)
                    sign = -1.0;
                else
                    sign = 1.0;
                
                A_matrix->val[i*2  ][j] = sign*pdcs.b[0][j];
                A_matrix->val[i*2+1][j] = sign*pdcs.b[1][j];
            }
            L_matrix->val[i*2    ][0]  = -pdcs.Fo;
            L_matrix->val[i*2 + 1][0]  = -pdcs.Go;
            
        }
        
        GMA_double_Tran(A_matrix,AT_matrix);
        GMA_double_mul(AT_matrix,A_matrix,ATA_matrix);
        GMA_double_inv(ATA_matrix,ATAI_matrix);
        GMA_double_mul(AT_matrix,L_matrix,ATL_matrix);
        GMA_double_mul(ATAI_matrix,ATL_matrix,X_matrix);
        
        GMA_double_mul(A_matrix,X_matrix,AX_matrix);
        GMA_double_sub(AX_matrix,L_matrix,V_matrix);
        
        double var_sum = 0;
        double max_V = -1000;
        for(int i=0;i<numofpts*2;i++)
        {
            /*
            //if(i%2 == 0)
            {
                double yy = V_matrix->val[i][0];
                printf("iter %d\tid %d\t%f\n",iter,i,yy);
            }
             */
            var_sum += (V_matrix->val[i][0])*(V_matrix->val[i][0]);
            if(max_V < fabs(V_matrix->val[i][0]))
                max_V = fabs(V_matrix->val[i][0]);
        }
        
        double sigma = sqrt(var_sum/(numofpts*2 - unknowns));
        
        eo.m_Wl += X_matrix->val[0][0];
        eo.m_Pl += X_matrix->val[1][0];
        eo.m_Kl += X_matrix->val[2][0];
        eo.m_Xl += X_matrix->val[3][0];
        eo.m_Yl += X_matrix->val[4][0];
        eo.m_Zl += X_matrix->val[5][0];
        
        double maxX = -100000;
        for(int i = 0;i<unknowns;i++)
        {
            if(maxX < fabs(X_matrix->val[i][0]))
                maxX = fabs(X_matrix->val[i][0]);
        }
        
        if(maxX < max_correct)
            check_stop = true;
        
        printf("iter %d\t eo %f\t%f\t%f\t%f\t%f\t%f\tsigma %f\tmax_V %f\tmaxX %f\n",iter,eo.m_Xl,eo.m_Yl,eo.m_Zl,
               eo.m_Wl,eo.m_Pl,eo.m_Kl,sigma,max_V,maxX);
    }
    
    GMA_double_destroy(A_matrix);
    GMA_double_destroy(L_matrix);
    GMA_double_destroy(AT_matrix);
    GMA_double_destroy(ATA_matrix);
    GMA_double_destroy(ATAI_matrix);
    GMA_double_destroy(ATL_matrix);
    GMA_double_destroy(X_matrix);
    GMA_double_destroy(AX_matrix);
    GMA_double_destroy(V_matrix);
}

void PPAEstimatefromEO(vector<D2DPOINT> &IPs, vector<D3DPOINT> &GCPs, EO &eo, CAMERA_INFO &camera)
{
    int numofpts = IPs.size();
    int unknowns = 2;
    
    GMA_double *A_matrix = GMA_double_create(2*numofpts,unknowns);
    GMA_double *L_matrix = GMA_double_create(2*numofpts,1);
    GMA_double *AT_matrix = GMA_double_create(unknowns,2*numofpts);
    GMA_double *ATA_matrix = GMA_double_create(unknowns,unknowns);
    GMA_double *ATAI_matrix = GMA_double_create(unknowns,unknowns);
    GMA_double *ATL_matrix  = GMA_double_create(unknowns,1);
    GMA_double *X_matrix    = GMA_double_create(unknowns,1);
    GMA_double *AX_matrix    = GMA_double_create(numofpts*2,1);
    GMA_double *V_matrix     = GMA_double_create(numofpts*2,1);
    
    bool check_stop = false;
    int max_iteration = 100;
    double max_correct = 0.000001;
    int iter = 0;
    while(!check_stop && iter < max_iteration)
    {
        iter++;
        
        UPARAMS uparam;
        PDCS pdcs;
        
        for(int i=0;i<numofpts;i++)
        {
            uparam.f = camera.m_focalLength;
            uparam.k1 = 0;
            uparam.k2 = 0;
            uparam.k3 = 0;
            uparam.p1 = 0;
            uparam.p2 = 0;
            uparam.xp = camera.m_ppx;
            uparam.yp = camera.m_ppy;
            uparam.a1 = 0;
            uparam.a2 = 0;
            
            uparam.X  = eo.m_Xl;
            uparam.Y  = eo.m_Yl;
            uparam.Z  = eo.m_Zl;
            uparam.omega = eo.m_Wl;
            uparam.phi   = eo.m_Pl;
            uparam.kappa = eo.m_Kl;
            
            uparam.xa = IPs[i].m_X;
            uparam.ya = IPs[i].m_Y;
            uparam.XA = GCPs[i].m_X;
            uparam.YA = GCPs[i].m_Y;
            uparam.ZA = GCPs[i].m_Z;
            
            GetPDCs(uparam,pdcs);
            
            A_matrix->val[i*2  ][0] = pdcs.c[0][1];
            A_matrix->val[i*2+1][0] = pdcs.c[1][1];
            
            A_matrix->val[i*2  ][1] = pdcs.c[0][2];
            A_matrix->val[i*2+1][1] = pdcs.c[1][2];
            
            L_matrix->val[i*2    ][0]  = -pdcs.Fo;
            L_matrix->val[i*2 + 1][0]  = -pdcs.Go;
            
        }
        
        GMA_double_Tran(A_matrix,AT_matrix);
        GMA_double_mul(AT_matrix,A_matrix,ATA_matrix);
        GMA_double_inv(ATA_matrix,ATAI_matrix);
        GMA_double_mul(AT_matrix,L_matrix,ATL_matrix);
        GMA_double_mul(ATAI_matrix,ATL_matrix,X_matrix);
        
        GMA_double_mul(A_matrix,X_matrix,AX_matrix);
        GMA_double_sub(AX_matrix,L_matrix,V_matrix);
        
        double var_sum = 0;
        double max_V = -1000;
        for(int i=0;i<numofpts*2;i++)
        {
            /*
            //if(i%2 == 0)
            {
                double yy = V_matrix->val[i][0];
                printf("iter %d\tid %d\t%f\n",iter,i,yy);
            }
             */
            var_sum += (V_matrix->val[i][0])*(V_matrix->val[i][0]);
            if(max_V < fabs(V_matrix->val[i][0]))
                max_V = fabs(V_matrix->val[i][0]);
        }
        
        double sigma = sqrt(var_sum/(numofpts*2 - unknowns));
        
        camera.m_ppx += X_matrix->val[0][0];
        camera.m_ppy += X_matrix->val[1][0];
        
        double maxX = -100000;
        for(int i = 0;i<unknowns;i++)
        {
            if(maxX < fabs(X_matrix->val[i][0]))
                maxX = fabs(X_matrix->val[i][0]);
        }
        
        if(maxX < max_correct)
            check_stop = true;
        
        printf("iter %d\t ppa %f\t%f\tsigma %f\tmax_V %f\tmaxX %f\n",iter,camera.m_ppx,camera.m_ppy,sigma,max_V,maxX);
    }
    
    GMA_double_destroy(A_matrix);
    GMA_double_destroy(L_matrix);
    GMA_double_destroy(AT_matrix);
    GMA_double_destroy(ATA_matrix);
    GMA_double_destroy(ATAI_matrix);
    GMA_double_destroy(ATL_matrix);
    GMA_double_destroy(X_matrix);
    GMA_double_destroy(AX_matrix);
    GMA_double_destroy(V_matrix);
}

void RadialParamEstimate(vector<D2DPOINT> &IPs, vector<D3DPOINT> &GCPs, EO &eo, CAMERA_INFO &camera)
{
    int numofpts = IPs.size();
    int unknowns = 3;
    
    GMA_double *A_matrix = GMA_double_create(2*numofpts,unknowns);
    GMA_double *L_matrix = GMA_double_create(2*numofpts,1);
    GMA_double *AT_matrix = GMA_double_create(unknowns,2*numofpts);
    GMA_double *ATA_matrix = GMA_double_create(unknowns,unknowns);
    GMA_double *ATAI_matrix = GMA_double_create(unknowns,unknowns);
    GMA_double *ATL_matrix  = GMA_double_create(unknowns,1);
    GMA_double *X_matrix    = GMA_double_create(unknowns,1);
    GMA_double *AX_matrix    = GMA_double_create(numofpts*2,1);
    GMA_double *V_matrix     = GMA_double_create(numofpts*2,1);
    
    bool check_stop = false;
    int max_iteration = 100;
    double max_correct = 0.000001;
    int iter = 0;
    while(!check_stop && iter < max_iteration)
    {
        iter++;
        
        UPARAMS uparam;
        PDCS pdcs;
        
        for(int i=0;i<numofpts;i++)
        {
            uparam.f = camera.m_focalLength;
            uparam.k1 = camera.k1;
            uparam.k2 = camera.k2;
            uparam.k3 = camera.k3;
            uparam.p1 = camera.p1;
            uparam.p2 = camera.p2;
            uparam.xp = camera.m_ppx;
            uparam.yp = camera.m_ppy;
            uparam.a1 = camera.a1;
            uparam.a2 = camera.a2;
            
            uparam.X  = eo.m_Xl;
            uparam.Y  = eo.m_Yl;
            uparam.Z  = eo.m_Zl;
            uparam.omega = eo.m_Wl;
            uparam.phi   = eo.m_Pl;
            uparam.kappa = eo.m_Kl;
            
            uparam.xa = IPs[i].m_X;
            uparam.ya = IPs[i].m_Y;
            uparam.XA = GCPs[i].m_X;
            uparam.YA = GCPs[i].m_Y;
            uparam.ZA = GCPs[i].m_Z;
            
            GetPDCs(uparam,pdcs);
            
            for(int j=0;j<2;j++)
            {
                A_matrix->val[i*2 + j][0] = pdcs.c[j][3];
                A_matrix->val[i*2 + j][1] = pdcs.c[j][4];
                A_matrix->val[i*2 + j][2] = pdcs.c[j][5];
            }
            
            L_matrix->val[i*2    ][0]  = -pdcs.Fo;
            L_matrix->val[i*2 + 1][0]  = -pdcs.Go;
            
        }
        
        GMA_double_Tran(A_matrix,AT_matrix);
        GMA_double_mul(AT_matrix,A_matrix,ATA_matrix);
        GMA_double_inv(ATA_matrix,ATAI_matrix);
        GMA_double_mul(AT_matrix,L_matrix,ATL_matrix);
        GMA_double_mul(ATAI_matrix,ATL_matrix,X_matrix);
        
        GMA_double_mul(A_matrix,X_matrix,AX_matrix);
        GMA_double_sub(AX_matrix,L_matrix,V_matrix);
        
        double var_sum = 0;
        double max_V = -1000;
        for(int i=0;i<numofpts*2;i++)
        {
            /*
            //if(i%2 == 0)
            {
                double yy = V_matrix->val[i][0];
                printf("iter %d\tid %d\t%f\n",iter,i,yy);
            }
             */
            var_sum += (V_matrix->val[i][0])*(V_matrix->val[i][0]);
            if(max_V < fabs(V_matrix->val[i][0]))
                max_V = fabs(V_matrix->val[i][0]);
        }
        
        double sigma = sqrt(var_sum/(numofpts*2 - unknowns));
        
        camera.k1 += X_matrix->val[0][0];
        camera.k2 += X_matrix->val[1][0];
        camera.k3 += X_matrix->val[2][0];
        
        double maxX = -100000;
        for(int i = 0;i<unknowns;i++)
        {
            if(maxX < fabs(X_matrix->val[i][0]))
                maxX = fabs(X_matrix->val[i][0]);
        }
        
        if(maxX < max_correct)
            check_stop = true;
        
        printf("iter %d\t k123 %f\t%f\t%f\tsigma %f\tmax_V %f\tmaxX %f\n",iter,camera.k1,camera.k2,camera.k3,sigma,max_V,maxX);
    }
    
    GMA_double_destroy(A_matrix);
    GMA_double_destroy(L_matrix);
    GMA_double_destroy(AT_matrix);
    GMA_double_destroy(ATA_matrix);
    GMA_double_destroy(ATAI_matrix);
    GMA_double_destroy(ATL_matrix);
    GMA_double_destroy(X_matrix);
    GMA_double_destroy(AX_matrix);
    GMA_double_destroy(V_matrix);
}

void TangentialParamEstimate(vector<D2DPOINT> &IPs, vector<D3DPOINT> &GCPs, EO &eo, CAMERA_INFO &camera)
{
    int numofpts = IPs.size();
    int unknowns = 2;
    
    GMA_double *A_matrix = GMA_double_create(2*numofpts,unknowns);
    GMA_double *L_matrix = GMA_double_create(2*numofpts,1);
    GMA_double *AT_matrix = GMA_double_create(unknowns,2*numofpts);
    GMA_double *ATA_matrix = GMA_double_create(unknowns,unknowns);
    GMA_double *ATAI_matrix = GMA_double_create(unknowns,unknowns);
    GMA_double *ATL_matrix  = GMA_double_create(unknowns,1);
    GMA_double *X_matrix    = GMA_double_create(unknowns,1);
    GMA_double *AX_matrix    = GMA_double_create(numofpts*2,1);
    GMA_double *V_matrix     = GMA_double_create(numofpts*2,1);
    
    bool check_stop = false;
    int max_iteration = 100;
    double max_correct = 0.000001;
    int iter = 0;
    while(!check_stop && iter < max_iteration)
    {
        iter++;
        
        UPARAMS uparam;
        PDCS pdcs;
        
        for(int i=0;i<numofpts;i++)
        {
            uparam.f = camera.m_focalLength;
            uparam.k1 = camera.k1;
            uparam.k2 = camera.k2;
            uparam.k3 = camera.k3;
            uparam.p1 = camera.p1;
            uparam.p2 = camera.p2;
            uparam.xp = camera.m_ppx;
            uparam.yp = camera.m_ppy;
            uparam.a1 = camera.a1;
            uparam.a2 = camera.a2;
            
            uparam.X  = eo.m_Xl;
            uparam.Y  = eo.m_Yl;
            uparam.Z  = eo.m_Zl;
            uparam.omega = eo.m_Wl;
            uparam.phi   = eo.m_Pl;
            uparam.kappa = eo.m_Kl;
            
            uparam.xa = IPs[i].m_X;
            uparam.ya = IPs[i].m_Y;
            uparam.XA = GCPs[i].m_X;
            uparam.YA = GCPs[i].m_Y;
            uparam.ZA = GCPs[i].m_Z;
            
            GetPDCs(uparam,pdcs);
            
            for(int j=0;j<2;j++)
            {
                A_matrix->val[i*2 + j][0] = pdcs.c[j][6];
                A_matrix->val[i*2 + j][1] = pdcs.c[j][7];
            }
            
            L_matrix->val[i*2    ][0]  = -pdcs.Fo;
            L_matrix->val[i*2 + 1][0]  = -pdcs.Go;
            
        }
        
        GMA_double_Tran(A_matrix,AT_matrix);
        GMA_double_mul(AT_matrix,A_matrix,ATA_matrix);
        GMA_double_inv(ATA_matrix,ATAI_matrix);
        GMA_double_mul(AT_matrix,L_matrix,ATL_matrix);
        GMA_double_mul(ATAI_matrix,ATL_matrix,X_matrix);
        
        GMA_double_mul(A_matrix,X_matrix,AX_matrix);
        GMA_double_sub(AX_matrix,L_matrix,V_matrix);
        
        double var_sum = 0;
        double max_V = -1000;
        for(int i=0;i<numofpts*2;i++)
        {
            /*
            //if(i%2 == 0)
            {
                double yy = V_matrix->val[i][0];
                printf("iter %d\tid %d\t%f\n",iter,i,yy);
            }
             */
            var_sum += (V_matrix->val[i][0])*(V_matrix->val[i][0]);
            if(max_V < fabs(V_matrix->val[i][0]))
                max_V = fabs(V_matrix->val[i][0]);
        }
        
        double sigma = sqrt(var_sum/(numofpts*2 - unknowns));
        
        camera.p1 += X_matrix->val[0][0];
        camera.p2 += X_matrix->val[1][0];
        
        double maxX = -100000;
        for(int i = 0;i<unknowns;i++)
        {
            if(maxX < fabs(X_matrix->val[i][0]))
                maxX = fabs(X_matrix->val[i][0]);
        }
        
        if(maxX < max_correct)
            check_stop = true;
        
        printf("iter %d\t p12 %f\t%f\tsigma %f\tmax_V %f\tmaxX %f\n",iter,camera.p1,camera.p2,sigma,max_V,maxX);
    }
    
    GMA_double_destroy(A_matrix);
    GMA_double_destroy(L_matrix);
    GMA_double_destroy(AT_matrix);
    GMA_double_destroy(ATA_matrix);
    GMA_double_destroy(ATAI_matrix);
    GMA_double_destroy(ATL_matrix);
    GMA_double_destroy(X_matrix);
    GMA_double_destroy(AX_matrix);
    GMA_double_destroy(V_matrix);
}

void AffinityParamEstimate(vector<D2DPOINT> &IPs, vector<D3DPOINT> &GCPs, EO &eo, CAMERA_INFO &camera)
{
    int numofpts = IPs.size();
    int unknowns = 2;
    
    GMA_double *A_matrix = GMA_double_create(2*numofpts,unknowns);
    GMA_double *L_matrix = GMA_double_create(2*numofpts,1);
    GMA_double *AT_matrix = GMA_double_create(unknowns,2*numofpts);
    GMA_double *ATA_matrix = GMA_double_create(unknowns,unknowns);
    GMA_double *ATAI_matrix = GMA_double_create(unknowns,unknowns);
    GMA_double *ATL_matrix  = GMA_double_create(unknowns,1);
    GMA_double *X_matrix    = GMA_double_create(unknowns,1);
    GMA_double *AX_matrix    = GMA_double_create(numofpts*2,1);
    GMA_double *V_matrix     = GMA_double_create(numofpts*2,1);
    
    bool check_stop = false;
    int max_iteration = 100;
    double max_correct = 0.000001;
    int iter = 0;
    while(!check_stop && iter < max_iteration)
    {
        iter++;
        
        UPARAMS uparam;
        PDCS pdcs;
        
        for(int i=0;i<numofpts;i++)
        {
            uparam.f = camera.m_focalLength;
            uparam.k1 = camera.k1;
            uparam.k2 = camera.k2;
            uparam.k3 = camera.k3;
            uparam.p1 = camera.p1;
            uparam.p2 = camera.p2;
            uparam.xp = camera.m_ppx;
            uparam.yp = camera.m_ppy;
            uparam.a1 = camera.a1;
            uparam.a2 = camera.a2;
            
            uparam.X  = eo.m_Xl;
            uparam.Y  = eo.m_Yl;
            uparam.Z  = eo.m_Zl;
            uparam.omega = eo.m_Wl;
            uparam.phi   = eo.m_Pl;
            uparam.kappa = eo.m_Kl;
            
            uparam.xa = IPs[i].m_X;
            uparam.ya = IPs[i].m_Y;
            uparam.XA = GCPs[i].m_X;
            uparam.YA = GCPs[i].m_Y;
            uparam.ZA = GCPs[i].m_Z;
            
            GetPDCs(uparam,pdcs);
            
            for(int j=0;j<2;j++)
            {
                A_matrix->val[i*2 + j][0] = pdcs.c[j][8];
                A_matrix->val[i*2 + j][1] = pdcs.c[j][9];
            }
            
            L_matrix->val[i*2    ][0]  = -pdcs.Fo;
            L_matrix->val[i*2 + 1][0]  = -pdcs.Go;
            
        }
        
        GMA_double_Tran(A_matrix,AT_matrix);
        GMA_double_mul(AT_matrix,A_matrix,ATA_matrix);
        GMA_double_inv(ATA_matrix,ATAI_matrix);
        GMA_double_mul(AT_matrix,L_matrix,ATL_matrix);
        GMA_double_mul(ATAI_matrix,ATL_matrix,X_matrix);
        
        GMA_double_mul(A_matrix,X_matrix,AX_matrix);
        GMA_double_sub(AX_matrix,L_matrix,V_matrix);
        
        double var_sum = 0;
        double max_V = -1000;
        for(int i=0;i<numofpts*2;i++)
        {
            /*
            //if(i%2 == 0)
            {
                double yy = V_matrix->val[i][0];
                printf("iter %d\tid %d\t%f\n",iter,i,yy);
            }
             */
            var_sum += (V_matrix->val[i][0])*(V_matrix->val[i][0]);
            if(max_V < fabs(V_matrix->val[i][0]))
                max_V = fabs(V_matrix->val[i][0]);
        }
        
        double sigma = sqrt(var_sum/(numofpts*2 - unknowns));
        
        camera.a1 += X_matrix->val[0][0];
        camera.a2 += X_matrix->val[1][0];
        
        double maxX = -100000;
        for(int i = 0;i<unknowns;i++)
        {
            if(maxX < fabs(X_matrix->val[i][0]))
                maxX = fabs(X_matrix->val[i][0]);
        }
        
        if(maxX < max_correct)
            check_stop = true;
        
        printf("iter %d\t a12 %f\t%f\tsigma %f\tmax_V %f\tmaxX %f\n",iter,camera.a1,camera.a2,sigma,max_V,maxX);
    }
    
    GMA_double_destroy(A_matrix);
    GMA_double_destroy(L_matrix);
    GMA_double_destroy(AT_matrix);
    GMA_double_destroy(ATA_matrix);
    GMA_double_destroy(ATAI_matrix);
    GMA_double_destroy(ATL_matrix);
    GMA_double_destroy(X_matrix);
    GMA_double_destroy(AX_matrix);
    GMA_double_destroy(V_matrix);
}

void Make_b_J_K(D3DPOINT PP, D3DPOINT GP, D2DPOINT p, EO eo, CAMERA_INFO camera, CLB &B)
{
    //Rotation Matrix
    RM M = MakeRotationMatrix(eo.m_Wl*RadToDeg, eo.m_Pl*RadToDeg, eo.m_Kl*RadToDeg);
    
    double Omega = eo.m_Wl;
    double Phi = eo.m_Pl;
    double Kappa = eo.m_Kl;
    
    double q, r, s;
    
    double dX, dY, dZ;
    
    double fqq;
    double f = camera.m_focalLength;
    D2DPOINT PPA(camera.m_ppx,camera.m_ppy);
    
    //printf("rot %f\t%f\t%f\ncamera %f\t%f\t%f\n",Omega,Phi,Kappa,f,PPA.m_X,PPA.m_Y);
    //printf("GP %f\t%f\t%f\tIP %f\t%f\n",GP.m_X,GP.m_Y,GP.m_Z,p.m_X,p.m_Y);
    dX = GP.m_X - PP.m_X;
    dY = GP.m_Y - PP.m_Y;
    dZ = GP.m_Z - PP.m_Z;

    q = M.m31*dX + M.m32*dY + M.m33*dZ;
    r = M.m11*dX + M.m12*dY + M.m13*dZ;
    s = M.m21*dX + M.m22*dY + M.m23*dZ;
    
    //printf("dxyz %f\t%f\t%f\t qrs %f\t%f\t%f\n",dX,dY,dZ,q,r,s);
    
    fqq = f/(q*q);

    B.b11 = fqq*(r*(-M.m33*dY+M.m32*dZ) - q*(-M.m13*dY+M.m12*dZ));
    B.b12 = fqq*(r*(cos(Phi)*dX+sin(Omega)*sin(Phi)*dY-cos(Omega)*sin(Phi)*dZ)
        - q*(-sin(Phi)*cos(Kappa)*dX+sin(Omega)*cos(Phi)*cos(Kappa)*dY-cos(Omega)*cos(Phi)*cos(Kappa)*dZ));
    B.b13 = -fqq*q*(M.m21*dX+M.m22*dY+M.m23*dZ);
    B.b14 = fqq*(r*M.m31-q*M.m11);
    B.b15 = fqq*(r*M.m32-q*M.m12);
    B.b16 = fqq*(r*M.m33-q*M.m13);
    
    B.J = p.m_X - PPA.m_X + f*r/q;
        
    B.b21 = fqq*(s*(-M.m33*dY+M.m32*dZ) - q*(-M.m23*dY+M.m22*dZ));
    B.b22 = fqq*(s*(cos(Phi)*dX+sin(Omega)*sin(Phi)*dY-cos(Omega)*sin(Phi)*dZ)
        - q*(sin(Phi)*sin(Kappa)*dX-sin(Omega)*cos(Phi)*sin(Kappa)*dY+cos(Omega)*cos(Phi)*sin(Kappa)*dZ));
    B.b23 = fqq*q*(M.m11*dX+M.m12*dY+M.m13*dZ);
    B.b24 = fqq*(s*M.m31-q*M.m21);
    B.b25 = fqq*(s*M.m32-q*M.m22);
    B.b26 = fqq*(s*M.m33-q*M.m23);

    B.K = p.m_Y - PPA.m_Y + f*s/q;
    /*
    printf("CLB\n%f\t%f\t%f\t%f\t%f\t%f\n",B.b11,B.b12,B.b13,B.b14,B.b15,B.b16);
    printf("%f\t%f\t%f\t%f\t%f\t%f\n",B.b21,B.b22,B.b23,B.b24,B.b25,B.b26);
    printf("%f\t%f\n",B.J, B.K);
    */
};

bool CalibrationBundle1(vector<D2DPOINT> &IPs, vector<D3DPOINT> &GCPs, EO &eo, CAMERA_INFO &camera)
{
    bool check_conver = false;
    
    int numofpts = IPs.size();
    
    GetInitialPCfromDLT(IPs, GCPs, eo, camera);
    
    int iteration = 0;
    int max_iter = 50;
    double max_correct = 10000000;
    
    //camera.m_ppx = 0;
    //camera.m_ppy = 0;
    //camera.m_focalLength = 8800.0;
    
    GMA_double *J_matrix = GMA_double_create(numofpts*2,16);
    GMA_double *K_matrix = GMA_double_create(numofpts*2,1);
    GMA_double *JT_matrix = GMA_double_create(16,numofpts*2);
    GMA_double *JTJ_matrix = GMA_double_create(16,16);
    GMA_double *JTJI_matrix = GMA_double_create(16,16);
    GMA_double *JTK_matrix = GMA_double_create(16,1);
    GMA_double *X_matrix = GMA_double_create(16,1);
    
    GMA_double *JX_matrix = GMA_double_create(numofpts*2,1);
    GMA_double *V_matrix = GMA_double_create(numofpts*2,1);
    
    printf("input camera %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",camera.m_focalLength,camera.m_ppx,camera.m_ppy,
           camera.k1,camera.k2, camera.k3,camera.p1,camera.p2,camera.a1,camera.a2);
    
    while(iteration < max_iter && !check_conver)
    {
        iteration++;
        
        UPARAMS uparam;
        PDCS pdcs;
        
        for(int i = 0 ; i < numofpts ; i++)
        {
            uparam.f  = camera.m_focalLength;
            uparam.k1 = camera.k1;
            uparam.k2 = camera.k2;
            uparam.k3 = camera.k3;
            uparam.p1 = camera.p1;
            uparam.p2 = camera.p2;
            uparam.xp = camera.m_ppx;
            uparam.yp = camera.m_ppy;
            uparam.a1 = camera.a1;
            uparam.a2 = camera.a2;
            
            uparam.omega    = eo.m_Wl;
            uparam.phi      = eo.m_Pl;
            uparam.kappa    = eo.m_Kl;
            uparam.X        = eo.m_Xl;
            uparam.Y        = eo.m_Yl;
            uparam.Z        = eo.m_Zl;
            
            uparam.xa = IPs[i].m_X;
            uparam.ya = IPs[i].m_Y;
            uparam.XA = GCPs[i].m_X;
            uparam.YA = GCPs[i].m_Y;
            uparam.ZA = GCPs[i].m_Z;
            
            //printf("uparam %f\t%f\t%f\t%f\t%f\t%f\t%f\n",uparam.f,uparam.X,uparam.Y,uparam.Z,uparam.omega,uparam.phi,uparam.kappa);
            
            
            //printf("IP XA %f\t%f\t%f\t%f\t%f\t%f\t%f\n",uparam.xa,uparam.ya,uparam.xp,uparam.yp,uparam.XA,uparam.YA,uparam.ZA);
            GetPDCs(uparam,pdcs);
            
            for(int j=0;j<6;j++)
            {
                double sign;
                if(j>2)
                    sign = -1.0;
                else
                    sign = 1.0;
                
                J_matrix->val[i*2  ][j] = sign*pdcs.b[0][j];
                J_matrix->val[i*2+1][j] = sign*pdcs.b[1][j];
            }
            
            for(int j=0; j<2; j++)
            {
                J_matrix->val[i*2+j][6 + 0] = pdcs.c[j][0];
                J_matrix->val[i*2+j][6 + 1] = pdcs.c[j][1];
                J_matrix->val[i*2+j][6 + 2] = pdcs.c[j][2];
                J_matrix->val[i*2+j][6 + 3] = pdcs.c[j][3];
                J_matrix->val[i*2+j][6 + 4] = pdcs.c[j][4];
                J_matrix->val[i*2+j][6 + 5] = pdcs.c[j][5];
                J_matrix->val[i*2+j][6 + 6] = pdcs.c[j][6];
                J_matrix->val[i*2+j][6 + 7] = pdcs.c[j][7];
                J_matrix->val[i*2+j][6 + 8] = pdcs.c[j][8];
                J_matrix->val[i*2+j][6 + 9] = pdcs.c[j][9];
            }
            
            
                    
            K_matrix->val[i*2  ][0] = -pdcs.Fo;
            K_matrix->val[i*2+1][0] = -pdcs.Go;
            
            //if(i > 4)
            //    exit(1);
        }
        //GMA_double_printf(K_matrix);
        
        
        
        GMA_double_Tran(J_matrix,JT_matrix);
        GMA_double_mul(JT_matrix,J_matrix,JTJ_matrix);
        GMA_double_inv(JTJ_matrix,JTJI_matrix);
        GMA_double_mul(JT_matrix,K_matrix,JTK_matrix);
        GMA_double_mul(JTJI_matrix,JTK_matrix,X_matrix);
        GMA_double_sub(JX_matrix,K_matrix,V_matrix);
        
        
        double var_sum = 0;
        double max_V = -1000;
        for(int i=0;i<numofpts*2;i++)
        {
            var_sum += (V_matrix->val[i][0])*(V_matrix->val[i][0]);
            if(max_V < fabs(V_matrix->val[i][0]))
                max_V = fabs(V_matrix->val[i][0]);
        }
        
        double sigma = sqrt(var_sum/(numofpts*2 - 16));
        printf("adjust sigma %f\n", sigma);
        //GMA_double_printf(X_matrix);
        
        //exit(1);
        
        eo.m_Wl += X_matrix->val[0][0];
        eo.m_Pl += X_matrix->val[1][0];
        eo.m_Kl += X_matrix->val[2][0];
        eo.m_Xl += X_matrix->val[3][0];
        eo.m_Yl += X_matrix->val[4][0];
        eo.m_Zl += X_matrix->val[5][0];
        
        camera.m_focalLength  += X_matrix->val[6][0];
        camera.m_ppx += X_matrix->val[7][0];
        camera.m_ppy += X_matrix->val[8][0];
        camera.k1 += X_matrix->val[9][0];
        camera.k2 += X_matrix->val[10][0];
        camera.k3 += X_matrix->val[11][0];
        camera.p1 += X_matrix->val[12][0];
        camera.p2 += X_matrix->val[13][0];
        camera.a1 += X_matrix->val[14][0];
        camera.a2 += X_matrix->val[15][0];
        
        printf("iter %d\teo %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",iteration,eo.m_Wl,eo.m_Pl,eo.m_Kl,eo.m_Xl,eo.m_Yl,eo.m_Zl,sigma,max_V);
        printf("iter %d\tca %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",iteration,camera.m_focalLength,camera.m_ppx,camera.m_ppy,
               camera.k1,camera.k2, camera.k3,camera.p1,camera.p2,camera.a1,camera.a2,
               sigma,max_V);
        
        printf("camera %f\t%f\t%f\n",camera.m_focalLength,camera.m_ppx,camera.m_ppy);
        
        max_correct = fabs(X_matrix->val[0][0]);
        for(int i=1;i<16;i++)
        {
            double t = X_matrix->val[i][0];
            printf("X %f\n",t);
            if(max_correct < fabs(X_matrix->val[i][0]))
                max_correct = fabs(X_matrix->val[i][0]);
        }
        
        if(max_correct <= 0.0000001)
            check_conver = true;
            
        printf("iter %d\tMax_correct %f\n",iteration,max_correct);
        
        //exit(1);
    }
    
    GMA_double_destroy(J_matrix);
    GMA_double_destroy(K_matrix);
    
    GMA_double_destroy(JT_matrix);
    GMA_double_destroy(JTJ_matrix);
    GMA_double_destroy(JTJI_matrix);
    GMA_double_destroy(JTK_matrix);
    GMA_double_destroy(X_matrix);
    GMA_double_destroy(JX_matrix);
    GMA_double_destroy(V_matrix);
    
    printf("Convergence %d\n",check_conver);
    
    return check_conver;
}

bool CalibrationBundlewithoutFL(vector<D2DPOINT> &IPs, vector<D3DPOINT> &GCPs, EO &eo, CAMERA_INFO &camera)
{
    bool check_conver = false;
    
    int numofpts = IPs.size();
    
    GetInitialPCfromDLT(IPs, GCPs, eo, camera);
    
    int iteration = 0;
    int max_iter = 100;
    double max_correct = 10000000;
    
    //camera.m_ppx = 0;
    //camera.m_ppy = 0;
    //camera.m_focalLength = 8800.0;
    int unknowns = 15;
    GMA_double *J_matrix = GMA_double_create(numofpts*2,unknowns);
    GMA_double *K_matrix = GMA_double_create(numofpts*2,1);
    GMA_double *JT_matrix = GMA_double_create(unknowns,numofpts*2);
    GMA_double *JTJ_matrix = GMA_double_create(unknowns,unknowns);
    GMA_double *JTJI_matrix = GMA_double_create(unknowns,unknowns);
    GMA_double *JTK_matrix = GMA_double_create(unknowns,1);
    GMA_double *X_matrix = GMA_double_create(unknowns,1);
    
    GMA_double *JX_matrix = GMA_double_create(numofpts*2,1);
    GMA_double *V_matrix = GMA_double_create(numofpts*2,1);
    
    while(iteration < max_iter && !check_conver)
    {
        iteration++;
        
        UPARAMS uparam;
        PDCS pdcs;
        
        for(int i = 0 ; i < numofpts ; i++)
        {
            uparam.f  = camera.m_focalLength;
            uparam.k1 = camera.k1;
            uparam.k2 = camera.k2;
            uparam.k3 = camera.k3;
            uparam.p1 = camera.p1;
            uparam.p2 = camera.p2;
            uparam.xp = camera.m_ppx;
            uparam.yp = camera.m_ppy;
            uparam.a1 = camera.a1;
            uparam.a2 = camera.a2;
            
            uparam.omega    = eo.m_Wl;
            uparam.phi      = eo.m_Pl;
            uparam.kappa    = eo.m_Kl;
            uparam.X        = eo.m_Xl;
            uparam.Y        = eo.m_Yl;
            uparam.Z        = eo.m_Zl;
            
            uparam.xa = IPs[i].m_X;
            uparam.ya = IPs[i].m_Y;
            uparam.XA = GCPs[i].m_X;
            uparam.YA = GCPs[i].m_Y;
            uparam.ZA = GCPs[i].m_Z;
            
            //printf("uparam %f\t%f\t%f\t%f\t%f\t%f\t%f\n",uparam.f,uparam.X,uparam.Y,uparam.Z,uparam.omega,uparam.phi,uparam.kappa);
            
            
            //printf("IP XA %f\t%f\t%f\t%f\t%f\t%f\t%f\n",uparam.xa,uparam.ya,uparam.xp,uparam.yp,uparam.XA,uparam.YA,uparam.ZA);
            GetPDCs(uparam,pdcs);
            
            for(int j=0;j<6;j++)
            {
                double sign;
                if(j>2)
                    sign = -1.0;
                else
                    sign = 1.0;
                
                J_matrix->val[i*2  ][j] = sign*pdcs.b[0][j];
                J_matrix->val[i*2+1][j] = sign*pdcs.b[1][j];
            }
            
            for(int j=0; j<2; j++)
            {
                J_matrix->val[i*2+j][6 + 0] = pdcs.c[j][1];
                J_matrix->val[i*2+j][6 + 1] = pdcs.c[j][2];
                J_matrix->val[i*2+j][6 + 2] = pdcs.c[j][3];
                J_matrix->val[i*2+j][6 + 3] = pdcs.c[j][4];
                J_matrix->val[i*2+j][6 + 4] = pdcs.c[j][5];
                J_matrix->val[i*2+j][6 + 5] = pdcs.c[j][6];
                J_matrix->val[i*2+j][6 + 6] = pdcs.c[j][7];
                J_matrix->val[i*2+j][6 + 7] = pdcs.c[j][8];
                J_matrix->val[i*2+j][6 + 8] = pdcs.c[j][9];
            }
            
            
                    
            K_matrix->val[i*2  ][0] = -pdcs.Fo;
            K_matrix->val[i*2+1][0] = -pdcs.Go;
            
            //if(i > 4)
            //    exit(1);
        }
        //GMA_double_printf(K_matrix);
        
        
        
        GMA_double_Tran(J_matrix,JT_matrix);
        GMA_double_mul(JT_matrix,J_matrix,JTJ_matrix);
        GMA_double_inv(JTJ_matrix,JTJI_matrix);
        GMA_double_mul(JT_matrix,K_matrix,JTK_matrix);
        GMA_double_mul(JTJI_matrix,JTK_matrix,X_matrix);
        GMA_double_sub(JX_matrix,K_matrix,V_matrix);
        
        
        double var_sum = 0;
        double max_V = -1000;
        for(int i=0;i<numofpts*2;i++)
        {
            var_sum += (V_matrix->val[i][0])*(V_matrix->val[i][0]);
            if(max_V < fabs(V_matrix->val[i][0]))
                max_V = fabs(V_matrix->val[i][0]);
        }
        
        double sigma = sqrt(var_sum/(numofpts*2 - unknowns));
        //GMA_double_printf(X_matrix);
        
        //exit(1);
        
        eo.m_Wl += X_matrix->val[0][0];
        eo.m_Pl += X_matrix->val[1][0];
        eo.m_Kl += X_matrix->val[2][0];
        eo.m_Xl += X_matrix->val[3][0];
        eo.m_Yl += X_matrix->val[4][0];
        eo.m_Zl += X_matrix->val[5][0];
        
        camera.m_ppx += X_matrix->val[6][0];
        camera.m_ppy += X_matrix->val[7][0];
        camera.k1 += X_matrix->val[8][0];
        camera.k2 += X_matrix->val[9][0];
        camera.k3 += X_matrix->val[10][0];
        camera.p1 += X_matrix->val[11][0];
        camera.p2 += X_matrix->val[12][0];
        camera.a1 += X_matrix->val[13][0];
        camera.a2 += X_matrix->val[14][0];
        
        printf("iter %d\teo %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",iteration,eo.m_Wl,eo.m_Pl,eo.m_Kl,eo.m_Xl,eo.m_Yl,eo.m_Zl,sigma,max_V);
        printf("iter %d\tca %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",iteration,camera.m_focalLength,camera.m_ppx,camera.m_ppy,
               camera.k1,camera.k2, camera.k3,camera.p1,camera.p2,camera.a1,camera.a2,
               sigma,max_V);
                
        //printf("camera %f\t%f\t%f\n",camera.m_focalLength,camera.m_ppx,camera.m_ppy);
        
        max_correct = fabs(X_matrix->val[0][0]);
        for(int i=1;i<unknowns;i++)
        {
            if(max_correct < fabs(X_matrix->val[i][0]))
                max_correct = fabs(X_matrix->val[i][0]);
        }
        
        printf("iter %d\tMax_correct %f\n",iteration,max_correct);
        
        if( max_correct <= 0.0001)
            check_conver = true;
    }
    
    GMA_double_destroy(J_matrix);
    GMA_double_destroy(K_matrix);
    
    GMA_double_destroy(JT_matrix);
    GMA_double_destroy(JTJ_matrix);
    GMA_double_destroy(JTJI_matrix);
    GMA_double_destroy(JTK_matrix);
    GMA_double_destroy(X_matrix);
    GMA_double_destroy(JX_matrix);
    GMA_double_destroy(V_matrix);
    
    printf("Convergence %d\n",check_conver);
    return check_conver;
}

void CalibrationBundle2(vector<D2DPOINT> &IPs, vector<D3DPOINT> &GCPs, EO &eo, CAMERA_INFO &camera)
{
    int numofpts = IPs.size();
    
    GetInitialPCfromDLT(IPs, GCPs, eo, camera);
    
    int iteration = 0;
    int max_iter = 100;
    double max_correct = 10000000;
    
    //camera.m_ppx = 0;
    //camera.m_ppy = 0;
    //camera.m_focalLength = 8800.0;
    
    GMA_double *J_matrix = GMA_double_create(numofpts*2,6);
    GMA_double *K_matrix = GMA_double_create(numofpts*2,1);
    GMA_double *JT_matrix = GMA_double_create(6,numofpts*2);
    GMA_double *JTJ_matrix = GMA_double_create(6,6);
    GMA_double *JTJI_matrix = GMA_double_create(6,6);
    GMA_double *JTK_matrix = GMA_double_create(6,1);
    GMA_double *X_matrix = GMA_double_create(6,1);
    
    GMA_double *JX_matrix = GMA_double_create(numofpts*2,1);
    GMA_double *V_matrix = GMA_double_create(numofpts*2,1);
    
    while(iteration < max_iter && max_correct > 0.00001)
    {
        iteration++;
        
        UPARAMS uparam;
        PDCS pdcs;
        
        for(int i = 0 ; i < numofpts ; i++)
        {
            uparam.f = camera.m_focalLength;
            uparam.k1 = 0;
            uparam.k2 = 0;
            uparam.k3 = 0;
            uparam.p1 = 0;
            uparam.p2 = 0;
            uparam.xp = camera.m_ppx;
            uparam.yp = camera.m_ppy;
            uparam.a1 = 0;
            uparam.a2 = 0;
            
            uparam.X  = eo.m_Xl;
            uparam.Y  = eo.m_Yl;
            uparam.Z  = eo.m_Zl;
            uparam.omega = eo.m_Wl;
            uparam.phi   = eo.m_Pl;
            uparam.kappa = eo.m_Kl;
            
            uparam.xa = IPs[i].m_X;
            uparam.ya = IPs[i].m_Y;
            uparam.XA = GCPs[i].m_X;
            uparam.YA = GCPs[i].m_Y;
            uparam.ZA = GCPs[i].m_Z;
            
            //printf("uparam %f\t%f\t%f\t%f\t%f\t%f\t%f\n",uparam.f,uparam.X,uparam.Y,uparam.Z,uparam.omega,uparam.phi,uparam.kappa);
            
            
            //printf("IP XA %f\t%f\t%f\t%f\t%f\t%f\t%f\n",uparam.xa,uparam.ya,uparam.xp,uparam.yp,uparam.XA,uparam.YA,uparam.ZA);
            GetPDCs(uparam,pdcs);
            
            for(int j=0;j<6;j++)
            {
                double sign;
                if(j>2)
                    sign = -1.0;
                else
                    sign = 1.0;
                
                J_matrix->val[i*2  ][j] = sign*pdcs.b[0][j];
                J_matrix->val[i*2+1][j] = sign*pdcs.b[1][j];
            }
            
            K_matrix->val[i*2  ][0] = -pdcs.Fo;
            K_matrix->val[i*2+1][0] = -pdcs.Go;
            
            //if(i > 4)
            //    exit(1);
        }
        //GMA_double_printf(K_matrix);
        
        
        
        GMA_double_Tran(J_matrix,JT_matrix);
        GMA_double_mul(JT_matrix,J_matrix,JTJ_matrix);
        GMA_double_inv(JTJ_matrix,JTJI_matrix);
        GMA_double_mul(JT_matrix,K_matrix,JTK_matrix);
        GMA_double_mul(JTJI_matrix,JTK_matrix,X_matrix);
        GMA_double_sub(JX_matrix,K_matrix,V_matrix);
        
        
        double var_sum = 0;
        double max_V = -1000;
        for(int i=0;i<numofpts*2;i++)
        {
            var_sum += (V_matrix->val[i][0])*(V_matrix->val[i][0]);
            if(max_V < fabs(V_matrix->val[i][0]))
                max_V = fabs(V_matrix->val[i][0]);
        }
        
        double sigma = sqrt(var_sum/(numofpts*2 - 6));
        //GMA_double_printf(X_matrix);
        
        //exit(1);
        
        eo.m_Wl += X_matrix->val[0][0];
        eo.m_Pl += X_matrix->val[1][0];
        eo.m_Kl += X_matrix->val[2][0];
        eo.m_Xl += X_matrix->val[3][0];
        eo.m_Yl += X_matrix->val[4][0];
        eo.m_Zl += X_matrix->val[5][0];
        
        
        
        max_correct = fabs(X_matrix->val[0][0]);
        for(int i=1;i<6;i++)
        {
            if(max_correct < fabs(X_matrix->val[i][0]))
                max_correct = fabs(X_matrix->val[i][0]);
        }
        
        printf("iter %d\teo %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",iteration,eo.m_Wl,eo.m_Pl,eo.m_Kl,eo.m_Xl,eo.m_Yl,eo.m_Zl,sigma,max_V);
        
        printf("camera %f\t%f\t%f\t%f\n",camera.m_focalLength,camera.m_ppx,camera.m_ppy,max_correct);
        
    }
    
    GMA_double_destroy(J_matrix);
    GMA_double_destroy(K_matrix);
    
    GMA_double_destroy(JT_matrix);
    GMA_double_destroy(JTJ_matrix);
    GMA_double_destroy(JTJI_matrix);
    GMA_double_destroy(JTK_matrix);
    GMA_double_destroy(X_matrix);
    GMA_double_destroy(JX_matrix);
    GMA_double_destroy(V_matrix);
    
    
    eo.m_Wl *= RadToDeg;
    eo.m_Pl *= RadToDeg;
    eo.m_Kl *= RadToDeg;
}

void CalibrationBundle(vector<D2DPOINT> &IPs, vector<D3DPOINT> &GCPs, EO &eo, CAMERA_INFO &camera)
{
    int numofpts = IPs.size();
    
    GetInitialPCfromDLT(IPs, GCPs, eo, camera);
    
    int iteration = 0;
    int max_iter = 100;
    double max_correct = 10000000;
    
    //camera.m_ppx = 0;
    //camera.m_ppy = 0;
    //camera.m_focalLength = 8800.0;
    
    GMA_double *J_matrix = GMA_double_create(numofpts*2,6+3);
    GMA_double *K_matrix = GMA_double_create(numofpts*2,1);
    GMA_double *JT_matrix = GMA_double_create(6+3,numofpts*2);
    GMA_double *JTJ_matrix = GMA_double_create(6+3,6+3);
    GMA_double *JTJI_matrix = GMA_double_create(6+3,6+3);
    GMA_double *JTK_matrix = GMA_double_create(6+3,1);
    GMA_double *X_matrix = GMA_double_create(6+3,1);
    
    GMA_double *JX_matrix = GMA_double_create(numofpts*2,1);
    GMA_double *V_matrix = GMA_double_create(numofpts*2,1);
    
    while(iteration < max_iter && max_correct > 0.0001)
    {
        iteration++;
        
        UPARAMS uparam;
        PDCS pdcs;
        
        for(int i = 0 ; i < numofpts ; i++)
        {
            uparam.f = camera.m_focalLength;
            uparam.k1 = 0;
            uparam.k2 = 0;
            uparam.k3 = 0;
            uparam.p1 = 0;
            uparam.p2 = 0;
            uparam.xp = camera.m_ppx;
            uparam.yp = camera.m_ppy;
            uparam.a1 = 0;
            uparam.a2 = 0;
            
            uparam.X  = eo.m_Xl;
            uparam.Y  = eo.m_Yl;
            uparam.Z  = eo.m_Zl;
            uparam.omega = eo.m_Wl;
            uparam.phi   = eo.m_Pl;
            uparam.kappa = eo.m_Kl;
            
            uparam.xa = IPs[i].m_X;
            uparam.ya = IPs[i].m_Y;
            uparam.XA = GCPs[i].m_X;
            uparam.YA = GCPs[i].m_Y;
            uparam.ZA = GCPs[i].m_Z;
            
            //printf("uparam %f\t%f\t%f\t%f\t%f\t%f\t%f\n",uparam.f,uparam.X,uparam.Y,uparam.Z,uparam.omega,uparam.phi,uparam.kappa);
            
            
            //printf("IP XA %f\t%f\t%f\t%f\t%f\t%f\t%f\n",uparam.xa,uparam.ya,uparam.xp,uparam.yp,uparam.XA,uparam.YA,uparam.ZA);
            GetPDCs(uparam,pdcs);
            
            for(int j=0;j<6;j++)
            {
                double sign;
                if(j>2)
                    sign = -1.0;
                else
                    sign = 1.0;
                
                J_matrix->val[i*2  ][j] = sign*pdcs.b[0][j];
                J_matrix->val[i*2+1][j] = sign*pdcs.b[1][j];
            }
            
            J_matrix->val[i*2  ][6 + 0] = pdcs.c[0][0];
            J_matrix->val[i*2+1][6 + 0] = pdcs.c[1][0];

            J_matrix->val[i*2  ][6 + 1] = pdcs.c[0][1];
            J_matrix->val[i*2+1][6 + 1] = pdcs.c[1][1];

            J_matrix->val[i*2  ][6 + 2] = pdcs.c[0][2];
            J_matrix->val[i*2+1][6 + 2] = pdcs.c[1][2];
                    
            K_matrix->val[i*2  ][0] = -pdcs.Fo;
            K_matrix->val[i*2+1][0] = -pdcs.Go;
            
            //if(i > 4)
            //    exit(1);
        }
        //GMA_double_printf(K_matrix);
        
        
        
        GMA_double_Tran(J_matrix,JT_matrix);
        GMA_double_mul(JT_matrix,J_matrix,JTJ_matrix);
        GMA_double_inv(JTJ_matrix,JTJI_matrix);
        GMA_double_mul(JT_matrix,K_matrix,JTK_matrix);
        GMA_double_mul(JTJI_matrix,JTK_matrix,X_matrix);
        GMA_double_sub(JX_matrix,K_matrix,V_matrix);
        
        
        double var_sum = 0;
        double max_V = -1000;
        for(int i=0;i<numofpts*2;i++)
        {
            var_sum += (V_matrix->val[i][0])*(V_matrix->val[i][0]);
            if(max_V < fabs(V_matrix->val[i][0]))
                max_V = fabs(V_matrix->val[i][0]);
        }
        
        double sigma = sqrt(var_sum/(numofpts*2 - 6));
        //GMA_double_printf(X_matrix);
        
        //exit(1);
        
        eo.m_Wl += X_matrix->val[0][0];
        eo.m_Pl += X_matrix->val[1][0];
        eo.m_Kl += X_matrix->val[2][0];
        eo.m_Xl += X_matrix->val[3][0];
        eo.m_Yl += X_matrix->val[4][0];
        eo.m_Zl += X_matrix->val[5][0];
        
        printf("iter %d\teo %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",iteration,eo.m_Wl,eo.m_Pl,eo.m_Kl,eo.m_Xl,eo.m_Yl,eo.m_Zl,sigma,max_V);
        
        camera.m_focalLength += X_matrix->val[6][0];
        camera.m_ppx += X_matrix->val[7][0];
        camera.m_ppy += X_matrix->val[8][0];
        
        printf("camera %f\t%f\t%f\n",camera.m_focalLength,camera.m_ppx,camera.m_ppy);
        
        max_correct =  -100000;
        for(int i=0;i<6+3;i++)
        {
            if(max_correct < fabs(X_matrix->val[i][0]))
                max_correct = fabs(X_matrix->val[i][0]);
            double tt = fabs(X_matrix->val[i][0]);
            //printf("X %d\t%f\n",i,tt);
        }
        
        printf("iter %d\tMax_correct %f\n",iteration,max_correct);
        
    }
    //printf("delete 1\n");
    GMA_double_destroy(J_matrix);
    GMA_double_destroy(K_matrix);
    
    GMA_double_destroy(JT_matrix);
    GMA_double_destroy(JTJ_matrix);
    GMA_double_destroy(JTJI_matrix);
    GMA_double_destroy(JTK_matrix);
    GMA_double_destroy(X_matrix);
    GMA_double_destroy(JX_matrix);
    GMA_double_destroy(V_matrix);
    
    //printf("initialize 1\n");
    J_matrix = GMA_double_create(numofpts*2,7);
    K_matrix = GMA_double_create(numofpts*2,1);
    JT_matrix = GMA_double_create(7,numofpts*2);
    JTJ_matrix = GMA_double_create(7,7);
    JTJI_matrix = GMA_double_create(7,7);
    JTK_matrix = GMA_double_create(7,1);
    X_matrix = GMA_double_create(7,1);
    
    JX_matrix = GMA_double_create(numofpts*2,1);
    V_matrix = GMA_double_create(numofpts*2,1);
    
    
    max_correct = 10000000;
    iteration = 0;
    while(iteration < max_iter && max_correct > 0.00001)
    {
        
        iteration++;
        UPARAMS uparam;
        PDCS pdcs;
        
        for(int i = 0 ; i < numofpts ; i++)
        {
            uparam.f  = camera.m_focalLength;
            uparam.k1 = camera.k1;
            uparam.k2 = camera.k2;
            uparam.k3 = camera.k3;
            uparam.p1 = camera.p1;
            uparam.p2 = camera.p2;
            uparam.xp = camera.m_ppx;
            uparam.yp = camera.m_ppy;
            uparam.a1 = camera.a1;
            uparam.a2 = camera.a2;
            
            uparam.omega = eo.m_Wl;
            uparam.phi = eo.m_Pl;
            uparam.kappa = eo.m_Kl;
            uparam.X = eo.m_Xl;
            uparam.Y = eo.m_Yl;
            uparam.Z = eo.m_Zl;
            
            uparam.xa = IPs[i].m_X;
            uparam.ya = IPs[i].m_Y;
            uparam.XA = GCPs[i].m_X;
            uparam.YA = GCPs[i].m_Y;
            uparam.ZA = GCPs[i].m_Z;
            
            //printf("initialize \n");
            
            GetPDCs(uparam,pdcs);
            //printf("GetPDCs \n");
            for(int j=0;j<2;j++)
            {
                J_matrix->val[i*2 + j][0] = pdcs.c[j][3];
                J_matrix->val[i*2 + j][1] = pdcs.c[j][4];

                J_matrix->val[i*2 + j][2] = pdcs.c[j][5];
                J_matrix->val[i*2 + j][3] = pdcs.c[j][6];

                J_matrix->val[i*2 + j][4] = pdcs.c[j][7];
                J_matrix->val[i*2 + j][5] = pdcs.c[j][8];
                J_matrix->val[i*2 + j][6] = pdcs.c[j][9];
            }
            //printf("J_matrix\n");
            K_matrix->val[i*2  ][0] = -pdcs.Fo;
            K_matrix->val[i*2+1][0] = -pdcs.Go;
            //printf("K_matrix\n");
        }
        //printf("matrix computation\n");
        //exit(1);
        GMA_double_Tran(J_matrix,JT_matrix);
        GMA_double_mul(JT_matrix,J_matrix,JTJ_matrix);
        GMA_double_inv(JTJ_matrix,JTJI_matrix);
        GMA_double_mul(JT_matrix,K_matrix,JTK_matrix);
        GMA_double_mul(JTJI_matrix,JTK_matrix,X_matrix);
        GMA_double_sub(JX_matrix,K_matrix,V_matrix);
        
        
        double var_sum = 0;
        double max_V = -1000;
        for(int i=0;i<numofpts*2;i++)
        {
            var_sum += (V_matrix->val[i][0])*(V_matrix->val[i][0]);
            if(max_V < fabs(V_matrix->val[i][0]))
                max_V = fabs(V_matrix->val[i][0]);
        }
        
        double sigma = sqrt(var_sum/(numofpts*2 - 6));
        
        camera.k1 += X_matrix->val[0][0];
        camera.k2 += X_matrix->val[1][0];
        camera.k3 += X_matrix->val[2][0];
        camera.p1 += X_matrix->val[3][0];
        camera.p2 += X_matrix->val[4][0];
        camera.a1 += X_matrix->val[5][0];
        camera.a2 += X_matrix->val[6][0];
        
        max_correct = fabs(X_matrix->val[0][0]);
        for(int i=1;i<7;i++)
        {
            if(max_correct < fabs(X_matrix->val[i][0]))
                max_correct = fabs(X_matrix->val[i][0]);
        }
        
        //printf("iter %d\t sigma %f\t param %f\t%f\t%f\t%f\t%f\t%f\t%f\n",iteration,sigma,camera.k1,camera.k2,camera.k3,camera.p1,camera.p2,camera.a1,camera.a2);
         
         
    }
     
    GMA_double_destroy(J_matrix);
    GMA_double_destroy(K_matrix);
    
    GMA_double_destroy(JT_matrix);
    GMA_double_destroy(JTJ_matrix);
    GMA_double_destroy(JTJI_matrix);
    GMA_double_destroy(JTK_matrix);
    GMA_double_destroy(X_matrix);
    GMA_double_destroy(JX_matrix);
    GMA_double_destroy(V_matrix);
    
    eo.m_Wl *= RadToDeg;
    eo.m_Pl *= RadToDeg;
    eo.m_Kl *= RadToDeg;
}

void GetPDCs(UPARAMS param, PDCS &PD)
{
    double RR;
    double x_, y_;
    x_ = param.xa-param.xp;
    y_ = param.ya-param.yp;
    RR = (x_*x_)+(y_*y_);

    //printf("x_ y_ %f\t%f\t%f\t%f\t%f\t%f\n",x_,y_,param.xa,param.ya,param.xp,param.yp);
    
    //dx
    PD.a[0][0] = param.k1*(RR+2*x_*x_) + param.k2*(RR*RR+4*RR*x_*x_) + param.k3*(RR*RR*RR+6*RR*RR*x_*x_)
            + 6*param.p1*x_ + 2*param.p2*y_ + 1;
    PD.a[1][0] = 2*param.k1*x_*y_ + 4*param.k2*RR*x_*y_ + 6*param.k3*RR*RR*x_*y_
            + 2*param.p1*y_ + 2*param.p2*x_ + param.a1;
    //dy
    PD.a[0][1] = 2*param.k1*x_*y_ + 4*param.k2*RR*x_*y_ + 6*param.k3*RR*RR*x_*y_
            + 2*param.p1*y_ + 2*param.p2*x_;
    PD.a[1][1] = param.k1*(RR+2*y_*y_) + param.k2*(RR*RR+4*RR*y_*y_) + param.k3*(RR*RR*RR+6*RR*RR*y_*y_)
            + 6*param.p2*y_ + 2*param.p1*x_ + param.a2 +1;

    //printf("PC.a %f\t%f\t%f\t%f\n",PD.a[0][0],PD.a[1][0],PD.a[0][1],PD.a[1][1]);
    
    double Rmatrix[3][3];
    Rmatrix[0][0] = cos(param.phi)*cos(param.kappa);
    Rmatrix[0][1] = sin(param.omega)*sin(param.phi)*cos(param.kappa) + cos(param.omega)*sin(param.kappa);
    Rmatrix[0][2] = -cos(param.omega)*sin(param.phi)*cos(param.kappa) + sin(param.omega)*sin(param.kappa);
    
    Rmatrix[1][0] = -cos(param.phi)*sin(param.kappa);
    Rmatrix[1][1] = -sin(param.omega)*sin(param.phi)*sin(param.kappa) + cos(param.omega)*cos(param.kappa);
    Rmatrix[1][2] = cos(param.omega)*sin(param.phi)*sin(param.kappa) + sin(param.omega)*cos(param.kappa);
    
    Rmatrix[2][0] = sin(param.phi);
    Rmatrix[2][1] = -sin(param.omega)*cos(param.phi);
    Rmatrix[2][2] = cos(param.omega)*cos(param.phi);
    /*
    printf("Rmatrix\n %f\t%f\t%f\n %f\t%f\t%f\n %f\t%f\t%f\n",
           Rmatrix[0][0],Rmatrix[0][1],Rmatrix[0][2],
           Rmatrix[1][0],Rmatrix[1][1],Rmatrix[1][2],
           Rmatrix[2][0],Rmatrix[2][1],Rmatrix[2][2]);
     */
    double q, r, s;
    double dX, dY, dZ;
    double fqq;

    dX = param.XA - param.X;
    dY = param.YA - param.Y;
    dZ = param.ZA - param.Z;

    //printf("dXYZ %f\t%f\t%f\n", dX,dY,dZ);
    
    q = Rmatrix[2][0]*dX + Rmatrix[2][1]*dY + Rmatrix[2][2]*dZ;
    if(q == 0.0)
        q = 1.0e-99;
    r = Rmatrix[0][0]*dX + Rmatrix[0][1]*dY + Rmatrix[0][2]*dZ;
    s = Rmatrix[1][0]*dX + Rmatrix[1][1]*dY + Rmatrix[1][2]*dZ;
    
    fqq = param.f/(q*q);

    //printf("qrsf %f\t%f\t%f\t%e\t%f\n",q,r,s,fqq,param.f);
    
    //omega, phi, kappa, X, Y, Z, (dXA, dYA, dZA)
    PD.b[0][0] = fqq*(q*(-Rmatrix[0][2]*dY+Rmatrix[0][1]*dZ) - r*(-Rmatrix[2][2]*dY+Rmatrix[2][1]*dZ));
    PD.b[0][1] = fqq*(q*(-sin(param.phi)*cos(param.kappa)*dX+sin(param.omega)*cos(param.phi)*cos(param.kappa)*dY-cos(param.omega)*cos(param.phi)*cos(param.kappa)*dZ)
            - r*(cos(param.phi)*dX+sin(param.omega)*sin(param.phi)*dY-cos(param.omega)*sin(param.phi)*dZ));
    PD.b[0][2] = fqq*q*(Rmatrix[1][0]*dX+Rmatrix[1][1]*dY+Rmatrix[1][2]*dZ);
    PD.b[0][3] = fqq*(q*Rmatrix[0][0] - r*Rmatrix[2][0]);
    PD.b[0][4] = fqq*(q*Rmatrix[0][1] - r*Rmatrix[2][1]);
    PD.b[0][5] = fqq*(q*Rmatrix[0][2] - r*Rmatrix[2][2]);

    PD.b[1][0] = fqq*(q*(-Rmatrix[1][2]*dY+Rmatrix[1][1]*dZ) - s*(-Rmatrix[2][2]*dY+Rmatrix[2][1]*dZ));
    PD.b[1][1] = fqq*(q*(sin(param.phi)*sin(param.kappa)*dX-sin(param.omega)*cos(param.phi)*sin(param.kappa)*dY+cos(param.omega)*cos(param.phi)*sin(param.kappa)*dZ)
            - s*(cos(param.phi)*dX+sin(param.omega)*sin(param.phi)*dY-cos(param.omega)*sin(param.phi)*dZ));
    PD.b[1][2] = -fqq*q*(Rmatrix[0][0]*dX+Rmatrix[0][1]*dY+Rmatrix[0][2]*dZ);
    PD.b[1][3] = fqq*(q*Rmatrix[1][0] - s*Rmatrix[2][0]);
    PD.b[1][4] = fqq*(q*Rmatrix[1][1] - s*Rmatrix[2][1]);
    PD.b[1][5] = fqq*(q*Rmatrix[1][2] - s*Rmatrix[2][2]);
    
    double dxR, dxT;
    double dyR, dyT;
    double dR_R;
    dR_R = param.k1*RR + param.k2*RR*RR + param.k3*RR*RR*RR;
    dxR = x_*dR_R;
    dyR = y_*dR_R;
    dxT = param.p1*(RR + 2*x_*x_) + 2*param.p2*x_*y_;
    dyT = param.p2*(RR + 2*y_*y_) + 2*param.p1*x_*y_;

    //calibration data (f, xp, yp, k1, k2, k3, p1, p2)
    PD.c[0][0] = r/q;
    PD.c[1][0] = s/q;

    PD.c[0][1] = -dR_R - 2*param.k1*x_*x_ - 4*param.k2*RR*x_*x_ - 6*param.k3*RR*RR*x_*x_
                - 6*param.p1*x_ - 2*param.p2*y_ - 1;
    PD.c[1][1] = -2*param.k1*x_*y_ - 4*param.k2*RR*x_*y_ - 6*param.k3*RR*RR*x_*y_
                - 2*param.p1*y_ - 2*param.p2*x_ - param.a1;
    PD.c[0][2] = -2*param.k1*x_*y_ - 4*param.k2*RR*x_*y_ - 6*param.k3*RR*RR*x_*y_
                - 2*param.p1*y_ - 2*param.p2*x_;
    PD.c[1][2] = -dR_R - 2*param.k1*y_*y_ - 4*param.k2*RR*y_*y_ - 6*param.k3*RR*RR*y_*y_
                - 6*param.p2*y_ - 2*param.p1*x_ - param.a2 - 1;

    PD.c[0][3] = x_*RR;
    PD.c[1][3] = y_*RR;
    PD.c[0][4] = x_*RR*RR;
    PD.c[1][4] = y_*RR*RR;
    PD.c[0][5] = x_*RR*RR*RR;
    PD.c[1][5] = y_*RR*RR*RR;

    PD.c[0][6] = RR + 2*x_*x_;
    PD.c[1][6] = 2*x_*y_;
    PD.c[0][7] = PD.c[1][6];
    PD.c[1][7] = RR + 2*y_*y_;

    PD.c[0][8] = 0.0;
    PD.c[1][8] = x_;
    PD.c[0][9] = 0.0;
    PD.c[1][9] = y_;

    //initial approximation of F and G
    PD.Fo = x_ + dxR + dxT + param.f*r/q;
    PD.Go = y_ + dyR + dyT + param.a1*x_ + param.a2*y_ + param.f*s/q;
}
