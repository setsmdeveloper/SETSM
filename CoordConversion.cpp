//
//  CoordConversion.cpp
//  
//
//  Created by Myoung-Jong Noh on 3/30/20.
//

#include "CoordConversion.hpp"


D2DPOINT *wgs2ps(TransParam _param, int _numofpts, D2DPOINT *_wgs)
{
    int m_NumOfPts = _numofpts;
    
    if(_param.projection == 1)
    {
        bool m_bHemisphere = _param.bHemisphere;
        double a = _param.a;
        double e = _param.e;
        double phi_c = _param.phi_c;
        double lambda_0 = _param.lambda_0;
        
        int pm = _param.pm;
        double t_c = _param.t_c;
        double m_c = _param.m_c;
        D2DPOINT *m_sPS;
        
        if (m_NumOfPts > 0) {
            m_sPS = (D2DPOINT *) malloc(sizeof(D2DPOINT) * m_NumOfPts);
        } else {
            return nullptr;
        }
        
#pragma omp parallel for schedule(guided)
        for (int i = 0; i < m_NumOfPts; i++) {
            double lambda = (double) (_wgs[i].m_X * pm * DegToRad);
            double phi= (double) (_wgs[i].m_Y * pm * DegToRad);
                
            double t = tan(PI / 4.0 - phi / 2.0) / pow((1.0 - e * sin(phi)) / (1.0 + e * sin(phi)), e / 2.0);
            double rho = a * m_c * t / t_c;
            
            double m = cos(phi) / sqrt(1.0 - pow(e, 2) * pow(sin(phi), 2));
            m_sPS[i].m_X = (double) (pm * rho * sin(lambda - lambda_0));
            m_sPS[i].m_Y = (double) (-pm * rho * cos(lambda - lambda_0));
        }
            
        return m_sPS;
    }
    else
    {
        double sa = _param.sa;
        double sb = _param.sb;
        double e2 = _param.e2;
        double e2cuadrada = _param.e2cuadrada;
        double c = _param.c;
        int Huso = _param.utm_zone;
        
        D2DPOINT *m_sPS;
        if (m_NumOfPts > 0) {
            m_sPS = (D2DPOINT *) malloc(sizeof(D2DPOINT) * m_NumOfPts);
        } else {
            return nullptr;
        }
        
#pragma omp parallel for schedule(guided)
        for (int i = 0; i < m_NumOfPts; i++) {
            double Lat = _wgs[i].m_Y;
            double Lon = _wgs[i].m_X;
            double lon = Lon * DegToRad;
            double lat = Lat * DegToRad;
            
            int S = ( ( Huso * 6 ) - 183 );
            double deltaS = lon - ( (double)S*DegToRad) ;
            
            double a = cos(lat)*sin(deltaS);
            double epsilon = 0.5 * log( ( 1 +  a)/ ( 1 - a ) );
            double nu = atan( tan(lat)/cos(deltaS) ) - lat;
            double v = ( c /  sqrt( 1 + ( e2cuadrada* cos(lat)*cos(lat)) ) )* 0.9996;
            double ta = ( e2cuadrada/ 2.0 ) * (epsilon*epsilon)* ( cos(lat)*cos(lat) );
            double a1 = sin( 2* lat );
            double a2 = a1* ( cos(lat)*cos(lat) );
            double j2 = lat + ( a1 / 2.0 );
            double j4 = ( ( 3 * j2 ) + a2 ) / 4.0;
            double j6 = ( ( 5 * j4 ) + ( a2 * ( cos(lat)*cos(lat) )) ) / 3;
            double alfa = ( 3 / 4.0 ) * e2cuadrada;
            double beta = ( 5 / 3.0 ) * alfa * alfa;
            double gama = ( 35 / 27.0 ) * alfa * alfa * alfa;
            double Bm = 0.9996 * c * ( lat - alfa * j2 + beta * j4 - gama * j6 );
            double xx = epsilon * v * ( 1 + ( ta / 3.0 ) ) + 500000;
            double yy = nu * v * ( 1 + ta ) + Bm;
            
            if (yy < 0)
                yy = 9999999 + yy;
            
            m_sPS[i].m_X = xx;
            m_sPS[i].m_Y = yy;
        }
        
        return m_sPS;
    }
}

D2DPOINT wgs2ps_single(TransParam _param, D2DPOINT _wgs)
{
    if(_param.projection == 1)
    {
        bool m_bHemisphere = _param.bHemisphere;
        double a = _param.a;
        double e = _param.e;
        double phi_c = _param.phi_c;
        double lambda_0 = _param.lambda_0;
        int pm = _param.pm;
        double t_c = _param.t_c;
        double m_c = _param.m_c;
        
        D2DPOINT m_sWGS = _wgs;
        D2DPOINT m_sPS;
        
        {
            m_sWGS.m_X = m_sWGS.m_X * pm * DegToRad;
            m_sWGS.m_Y = m_sWGS.m_Y * pm * DegToRad;
            double lambda = m_sWGS.m_X;
            double phi = m_sWGS.m_Y;
            
            double t = tan(PI / 4.0 - phi / 2.0) / pow((1.0 - e * sin(phi)) / (1.0 + e * sin(phi)), e / 2.0);
            double rho = a * m_c * t / t_c;
            
            double m = cos(phi) / sqrt(1.0 - pow(e, 2) * pow(sin(phi), 2));
            m_sPS.m_X = pm * rho * sin(lambda - lambda_0);
            m_sPS.m_Y = -pm * rho * cos(lambda - lambda_0);
        }
        
        return m_sPS;
    }
    else
    {
        double sa = _param.sa;
        double sb = _param.sb;
        double e2 = _param.e2;
        double e2cuadrada = _param.e2cuadrada;
        double c = _param.c;
        int Huso = _param.utm_zone;
        
        D2DPOINT m_sPS;
        
        double Lat = _wgs.m_Y;
        double Lon = _wgs.m_X;
        double lon = Lon * DegToRad;
        double lat = Lat * DegToRad;
        
        int S = ( ( Huso * 6 ) - 183 );
        double deltaS = lon - ( (double)S*DegToRad) ;
        
        double a = cos(lat)*sin(deltaS);
        double epsilon = 0.5 * log( ( 1 +  a)/ ( 1 - a ) );
        double nu = atan( tan(lat)/cos(deltaS) ) - lat;
        double v = ( c /  sqrt( 1 + ( e2cuadrada* cos(lat)*cos(lat)) ) )* 0.9996;
        double ta = ( e2cuadrada/ 2.0 ) * (epsilon*epsilon)* ( cos(lat)*cos(lat) );
        double a1 = sin( 2* lat );
        double a2 = a1* ( cos(lat)*cos(lat) );
        double j2 = lat + ( a1 / 2.0 );
        double j4 = ( ( 3 * j2 ) + a2 ) / 4.0;
        double j6 = ( ( 5 * j4 ) + ( a2 * ( cos(lat)*cos(lat) )) ) / 3;
        double alfa = ( 3 / 4.0 ) * e2cuadrada;
        double beta = ( 5 / 3.0 ) * alfa * alfa;
        double gama = ( 35 / 27.0 ) * alfa * alfa * alfa;
        double Bm = 0.9996 * c * ( lat - alfa * j2 + beta * j4 - gama * j6 );
        double xx = epsilon * v * ( 1 + ( ta / 3.0 ) ) + 500000;
        double yy = nu * v * ( 1 + ta ) + Bm;
        
        
        if (yy < 0)
            yy = 9999999 + yy;
        
        m_sPS.m_X = xx;
        m_sPS.m_Y = yy;
        
        return m_sPS;
    }
    
}

D3DPOINT wgs2ps_single_3D(TransParam _param, D3DPOINT _wgs)
{
    if(_param.projection == 1)
    {
        bool m_bHemisphere = _param.bHemisphere;
        double a = _param.a;
        double e = _param.e;
        double phi_c = _param.phi_c;
        double lambda_0 = _param.lambda_0;
        int pm = _param.pm;
        double t_c = _param.t_c;
        double m_c = _param.m_c;
        
        D2DPOINT m_sWGS;
        m_sWGS.m_X = _wgs.m_X;
        m_sWGS.m_Y = _wgs.m_Y;
        D3DPOINT m_sPS;
        
        {
            m_sWGS.m_X = m_sWGS.m_X * pm * DegToRad;
            m_sWGS.m_Y = m_sWGS.m_Y * pm * DegToRad;
            double lambda = m_sWGS.m_X;
            double phi = m_sWGS.m_Y;
            
            double t = tan(PI / 4.0 - phi / 2.0) / pow((1.0 - e * sin(phi)) / (1.0 + e * sin(phi)), e / 2.0);
            double rho = a * m_c * t / t_c;
            
            double m = cos(phi) / sqrt(1.0 - pow(e, 2) * pow(sin(phi), 2));
            m_sPS.m_X = pm * rho * sin(lambda - lambda_0);
            m_sPS.m_Y = -pm * rho * cos(lambda - lambda_0);
            m_sPS.m_Z = _wgs.m_Z;
        }
        
        return m_sPS;
    }
    else
    {
        double sa = _param.sa;
        double sb = _param.sb;
        double e2 = _param.e2;
        double e2cuadrada = _param.e2cuadrada;
        double c = _param.c;
        int Huso = _param.utm_zone;
        
        D3DPOINT m_sPS;
        
        double Lat = _wgs.m_Y;
        double Lon = _wgs.m_X;
        double lon = Lon * DegToRad;
        double lat = Lat * DegToRad;
        
        int S = ( ( Huso * 6 ) - 183 );
        double deltaS = lon - ( (double)S*DegToRad) ;
        
        double a = cos(lat)*sin(deltaS);
        double epsilon = 0.5 * log( ( 1 +  a)/ ( 1 - a ) );
        double nu = atan( tan(lat)/cos(deltaS) ) - lat;
        double v = ( c /  sqrt( 1 + ( e2cuadrada* cos(lat)*cos(lat)) ) )* 0.9996;
        double ta = ( e2cuadrada/ 2.0 ) * (epsilon*epsilon)* ( cos(lat)*cos(lat) );
        double a1 = sin( 2* lat );
        double a2 = a1* ( cos(lat)*cos(lat) );
        double j2 = lat + ( a1 / 2.0 );
        double j4 = ( ( 3 * j2 ) + a2 ) / 4.0;
        double j6 = ( ( 5 * j4 ) + ( a2 * ( cos(lat)*cos(lat) )) ) / 3;
        double alfa = ( 3 / 4.0 ) * e2cuadrada;
        double beta = ( 5 / 3.0 ) * alfa * alfa;
        double gama = ( 35 / 27.0 ) * alfa * alfa * alfa;
        double Bm = 0.9996 * c * ( lat - alfa * j2 + beta * j4 - gama * j6 );
        double xx = epsilon * v * ( 1 + ( ta / 3.0 ) ) + 500000;
        double yy = nu * v * ( 1 + ta ) + Bm;
        
        
        if (yy < 0)
            yy = 9999999 + yy;
        
        m_sPS.m_X = xx;
        m_sPS.m_Y = yy;
        m_sPS.m_Z = _wgs.m_Z;
        
        return m_sPS;
    }
}


D3DPOINT *wgs2ps_3D(TransParam _param, int _numofpts, D3DPOINT *_wgs)
{
    int m_NumOfPts = _numofpts;
    
    if(_param.projection == 1)
    {
        
        bool m_bHemisphere = _param.bHemisphere;
        double a = _param.a;
        double e = _param.e;
        double phi_c = _param.phi_c;
        double lambda_0 = _param.lambda_0;
        
        int pm = _param.pm;
        double t_c = _param.t_c;
        double m_c = _param.m_c;
        D3DPOINT *m_sPS;
        
        if (m_NumOfPts > 0) {
            m_sPS = (D3DPOINT *) malloc(sizeof(D3DPOINT) * m_NumOfPts);
        } else {
            return nullptr;
        }
        
#pragma omp parallel for schedule(guided)
        for (int i = 0; i < m_NumOfPts; i++) {
          double lambda = _wgs[i].m_X * pm * DegToRad;
          double phi = _wgs[i].m_Y * pm * DegToRad;
            
            double t = tan(PI / 4.0 - phi / 2.0) / pow((1.0 - e * sin(phi)) / (1.0 + e * sin(phi)), e / 2.0);
            double rho = a * m_c * t / t_c;
            
            double m = cos(phi) / sqrt(1.0 - pow(e, 2) * pow(sin(phi), 2));
            m_sPS[i].m_X = pm * rho * sin(lambda - lambda_0);
            m_sPS[i].m_Y = -pm * rho * cos(lambda - lambda_0);
            m_sPS[i].m_Z = _wgs[i].m_Z;
        }
            
        return m_sPS;
    }
    else
    {
        double sa = _param.sa;
        double sb = _param.sb;
        double e2 = _param.e2;
        double e2cuadrada = _param.e2cuadrada;
        double c = _param.c;
        int Huso = _param.utm_zone;
        
        D3DPOINT *m_sPS;
        
        if (m_NumOfPts > 0) {
            m_sPS = (D3DPOINT *) malloc(sizeof(D3DPOINT) * m_NumOfPts);
        } else {
            return nullptr;
        }
        
#pragma omp parallel for schedule(guided)
        for (int i = 0; i < m_NumOfPts; i++) {
            double Lat = _wgs[i].m_Y;
            double Lon = _wgs[i].m_X;
            double lon = Lon * DegToRad;
            double lat = Lat * DegToRad;
            
            int S = ( ( Huso * 6 ) - 183 );
            double deltaS = lon - ( (double)S*DegToRad) ;
            
            double a = cos(lat)*sin(deltaS);
            double epsilon = 0.5 * log( ( 1 +  a)/ ( 1 - a ) );
            double nu = atan( tan(lat)/cos(deltaS) ) - lat;
            double v = ( c /  sqrt( 1 + ( e2cuadrada* cos(lat)*cos(lat)) ) )* 0.9996;
            double ta = ( e2cuadrada/ 2.0 ) * (epsilon*epsilon)* ( cos(lat)*cos(lat) );
            double a1 = sin( 2* lat );
            double a2 = a1* ( cos(lat)*cos(lat) );
            double j2 = lat + ( a1 / 2.0 );
            double j4 = ( ( 3 * j2 ) + a2 ) / 4.0;
            double j6 = ( ( 5 * j4 ) + ( a2 * ( cos(lat)*cos(lat) )) ) / 3;
            double alfa = ( 3 / 4.0 ) * e2cuadrada;
            double beta = ( 5 / 3.0 ) * alfa * alfa;
            double gama = ( 35 / 27.0 ) * alfa * alfa * alfa;
            double Bm = 0.9996 * c * ( lat - alfa * j2 + beta * j4 - gama * j6 );
            double xx = epsilon * v * ( 1 + ( ta / 3.0 ) ) + 500000;
            double yy = nu * v * ( 1 + ta ) + Bm;
            
            
            if (yy < 0)
                yy = 9999999 + yy;
            
            m_sPS[i].m_X = xx;
            m_sPS[i].m_Y = yy;
            m_sPS[i].m_Z = _wgs[i].m_Z;
        }
        
        return m_sPS;
        
    }
    
}

D2DPOINT *ps2wgs(const TransParam _param,const long int _numofpts,const D2DPOINT *_ps)
{
    const long int m_NumOfPts = _numofpts;
   
    if(m_NumOfPts > 0)
    {
        if(_param.projection == 1)
        {
            const bool m_bHemisphere = _param.bHemisphere;
            const double a = _param.a;
            const double e = _param.e;
            const double phi_c = _param.phi_c;
            const double lambda_0 = _param.lambda_0;

            const int pm = _param.pm;
            const double t_c = _param.t_c;
            const double m_c = _param.m_c;

            D2DPOINT *m_sWGS = (D2DPOINT *) malloc(sizeof(D2DPOINT) * m_NumOfPts);

            const double e2 = e * e;
            const double e4 = e2 * e2;
            const double e6 = e2 * e4;
            const double e8 = e4 * e4;

#pragma omp parallel for schedule(guided)
            for (int i = 0; i < m_NumOfPts; i++)
            {
                const double x = _ps[i].m_X * pm;
                const double y = _ps[i].m_Y * pm;

                const double rho = sqrt(pow(x, 2) + pow(y, 2));
                const double t = rho * t_c / (a * m_c);

                const double chi = PI / 2 - 2 * atan(t);
                double phi = chi + (e2 / 2 + 5 * e4 / 24 + e6 / 12 + 13 * e8 / 360) * sin(2 * chi) + (7 * e4 / 48 + 29 * e6 / 240 + 811 * e8 / 11520) * sin(4 * chi) +
                    (7 * e6 / 120 + 81 * e8 / 1120) * sin(6 * chi) + (4279 * e8 / 161280) * sin(8 * chi);

                double lambda = lambda_0 + atan2(x, -y);
                phi = pm * phi;
                lambda = pm * lambda;

                if (lambda > PI) {
                    lambda = lambda - 2 * PI;
                } else if (lambda < -PI) {
                    lambda = lambda + 2 * PI;
                }

                m_sWGS[i].m_Y = RadToDeg * phi;
                m_sWGS[i].m_X = RadToDeg * lambda;
            }

            return m_sWGS;
        }
        else
        {
            const int hemis = _param.pm; //1 = north, -1 = south
            const double sa = _param.sa;
            const double sb = _param.sb;
            const double e2 = _param.e2;
            const double e2cuadrada = _param.e2cuadrada;
            const double c = _param.c;

            D2DPOINT *m_sWGS = (D2DPOINT *) malloc(sizeof(D2DPOINT) * m_NumOfPts);

#pragma omp parallel for schedule(guided)
            for (int i = 0; i < m_NumOfPts; i++)
            {
                const double x = _ps[i].m_X;
                const double y = _ps[i].m_Y;

                const double X = x - 500000;
                double Y = y;
                if (hemis < 0)
                    Y = Y - 10000000;

                const int S = ( ( _param.utm_zone* 6 ) - 183 );
                const double lat =  Y / ( 6366197.724 * 0.9996 );
                const double v = ( c / ( sqrt( 1.0 + ( e2cuadrada * ( cos(lat)*cos(lat) ) ) ) ) ) * 0.9996;
                const double a = X / v;


                const double a1 = sin( 2* lat );
                const double a2 = a1* ( cos(lat)*cos(lat) );
                const double j2 = lat + ( a1 / 2.0 );
                const double j4 = ( ( 3 * j2 ) + a2 ) / 4.0;
                const double j6 = ( ( 5 * j4 ) + ( a2 * ( cos(lat)*cos(lat) )) ) / 3;
                const double alfa = ( 3 / 4.0 ) * e2cuadrada;
                const double beta = ( 5 / 3.0 ) * alfa * alfa;
                const double gama = ( 35 / 27.0 ) * alfa * alfa * alfa;
                const double Bm = 0.9996 * c * ( lat - alfa * j2 + beta * j4 - gama * j6 );
                const double b = ( Y - Bm ) / v;
                const double Epsi = ( ( e2cuadrada * a *a ) / 2.0 ) * ( cos(lat)*cos(lat) );
                const double Eps = a * ( 1 - ( Epsi / 3.0 ) );
                const double nab = ( b * ( 1 - Epsi ) ) + lat;
                const double senoheps = ( exp(Eps) - exp(-Eps) ) / 2;
                const double Delt = atan(senoheps / (cos(nab) ) );
                const double TaO = atan(cos(Delt) * tan(nab));
                const double longitude = (Delt * (180/PI) ) + S;

                const double latitude = ( lat + ( 1 + e2cuadrada * (cos(lat)*cos(lat)) - ( 3/2.0 )* e2cuadrada * sin(lat) * cos(lat) * ( TaO - lat ) )* ( TaO - lat ) ) * (180/PI);

                m_sWGS[i].m_Y = latitude;
                m_sWGS[i].m_X = longitude;
            }
            return m_sWGS;
        }
    }
    else
    {
        return NULL;
    }
}

D2DPOINT ps2wgs_single(const TransParam _param, const D2DPOINT _ps)
{
    if(_param.projection == 1)
    {
        bool m_bHemisphere = _param.bHemisphere;
        double a = _param.a;
        double e = _param.e;
        double phi_c = _param.phi_c;
        double lambda_0 = _param.lambda_0;
        int pm = _param.pm;
        double t_c = _param.t_c;
        double m_c = _param.m_c;
        
        D2DPOINT m_sWGS;
        D2DPOINT m_sPS = _ps;
        
        double e2 = e * e;
        double e4 = e2 * e2;
        double e6 = e2 * e4;
        double e8 = e4 * e4;
        
        m_sPS.m_X = m_sPS.m_X * pm;
        m_sPS.m_Y = m_sPS.m_Y * pm;
        
        double rho = sqrt(pow(m_sPS.m_X, 2) + pow(m_sPS.m_Y, 2));
        double t = rho * t_c / (a * m_c);
        
        double chi = PI / 2 - 2 * atan(t);
        double phi = chi + (e2 / 2 + 5 * e4 / 24 + e6 / 12 + 13 * e8 / 360) * sin(2 * chi) + (7 * e4 / 48 + 29 * e6 / 240 + 811 * e8 / 11520) * sin(4 * chi) +
            (7 * e6 / 120 + 81 * e8 / 1120) * sin(6 * chi) + (4279 * e8 / 161280) * sin(8 * chi);
        
        double lambda = lambda_0 + atan2(m_sPS.m_X, -m_sPS.m_Y);
        phi = pm * phi;
        lambda = pm * lambda;
        
        if (lambda > PI) {
            lambda = lambda - 2 * PI;
        } else if (lambda < -PI) {
            lambda = lambda + 2 * PI;
        }
        
        m_sWGS.m_Y = RadToDeg * phi;
        m_sWGS.m_X = RadToDeg * lambda;
        
        return m_sWGS;
    }
    else
    {
        int hemis = _param.pm; //1 = north, -1 = south
//        printf("hemis %d\n",hemis);
        double sa = _param.sa;
        double sb = _param.sb;
        double e2 = _param.e2;
        double e2cuadrada = _param.e2cuadrada;
        double c = _param.c;
        
        D2DPOINT m_sWGS;
        D2DPOINT m_sPS = _ps;
        
        double x = m_sPS.m_X;
        double y = m_sPS.m_Y;
        
        double X = x - 500000;
        double Y = y;
        if (hemis < 0)
            Y = Y - 10000000;
        
        int S = ( ( _param.utm_zone* 6 ) - 183 );
        double lat =  Y / ( 6366197.724 * 0.9996 );
        double v = ( c / ( sqrt( 1.0 + ( e2cuadrada * ( cos(lat)*cos(lat) ) ) ) ) ) * 0.9996;
        double a = X / v;
        
        
        double a1 = sin( 2* lat );
        double a2 = a1* ( cos(lat)*cos(lat) );
        double j2 = lat + ( a1 / 2.0 );
        double j4 = ( ( 3 * j2 ) + a2 ) / 4.0;
        double j6 = ( ( 5 * j4 ) + ( a2 * ( cos(lat)*cos(lat) )) ) / 3;
        double alfa = ( 3 / 4.0 ) * e2cuadrada;
        double beta = ( 5 / 3.0 ) * alfa * alfa;
        double gama = ( 35 / 27.0 ) * alfa * alfa * alfa;
        double Bm = 0.9996 * c * ( lat - alfa * j2 + beta * j4 - gama * j6 );
        double b = ( Y - Bm ) / v;
        double Epsi = ( ( e2cuadrada * a *a ) / 2.0 ) * ( cos(lat)*cos(lat) );
        double Eps = a * ( 1 - ( Epsi / 3.0 ) );
        double nab = ( b * ( 1 - Epsi ) ) + lat;
        double senoheps = ( exp(Eps) - exp(-Eps) ) / 2;
        double Delt = atan(senoheps / (cos(nab) ) );
        double TaO = atan(cos(Delt) * tan(nab));
        double longitude = (Delt * (180/PI) ) + S;
        
        double latitude = ( lat + ( 1 + e2cuadrada * (cos(lat)*cos(lat)) - ( 3/2.0 )* e2cuadrada * sin(lat) * cos(lat) * ( TaO - lat ) )* ( TaO - lat ) ) * (180/PI);
        
        m_sWGS.m_Y = latitude;
        m_sWGS.m_X = longitude;
        
        return m_sWGS;
    }
    
}

D3DPOINT *ps2wgs_3D(TransParam _param, int _numofpts, D3DPOINT *_ps)
{
    int m_NumOfPts = _numofpts;
    if(_param.projection == 1)
    {
        bool m_bHemisphere = _param.bHemisphere;
        double a = _param.a;
        double e = _param.e;
        double phi_c = _param.phi_c;
        double lambda_0 = _param.lambda_0;
        
        int pm = _param.pm;
        double t_c = _param.t_c;
        double m_c = _param.m_c;
            
        D3DPOINT *m_sWGS;
        
        if (m_NumOfPts > 0) {
            m_sWGS = (D3DPOINT *) malloc(sizeof(D3DPOINT) * m_NumOfPts);
        } else {
            return nullptr;
        }
        
        double e2 = e * e;
        double e4 = e2 * e2;
        double e6 = e2 * e4;
        double e8 = e4 * e4;
        
#pragma omp parallel for schedule(guided)
        for (int i = 0; i < m_NumOfPts; i++) {
          double x = _ps[i].m_X * pm;
          double y = _ps[i].m_Y * pm;
            
            double rho = sqrt(pow(x, 2) + pow(y, 2));
            double t = rho * t_c / (a * m_c);
            
            double chi = PI / 2 - 2 * atan(t);
            double phi = chi + (e2 / 2 + 5 * e4 / 24 + e6 / 12 + 13 * e8 / 360) * sin(2 * chi) + (7 * e4 / 48 + 29 * e6 / 240 + 811 * e8 / 11520) * sin(4 * chi) +
                (7 * e6 / 120 + 81 * e8 / 1120) * sin(6 * chi) + (4279 * e8 / 161280) * sin(8 * chi);
            
            double lambda = lambda_0 + atan2(x, -y);
            phi = pm * phi;
            lambda = pm * lambda;
            if (lambda > PI) {
                lambda = lambda - 2 * PI;
            } else if (lambda < -PI) {
                lambda = lambda + 2 * PI;
            }
            
            m_sWGS[i].m_Y = RadToDeg * phi;
            m_sWGS[i].m_X = RadToDeg * lambda;
            m_sWGS[i].m_Z = _ps[i].m_Z;
        }
            
        return m_sWGS;
    }
    else
    {
        int hemis = _param.pm; //1 = north, -1 = south
        double sa = _param.sa;
        double sb = _param.sb;
        double e2 = _param.e2;
        double e2cuadrada = _param.e2cuadrada;
        double c = _param.c;
            
        D3DPOINT *m_sWGS;
        
        if (m_NumOfPts > 0) {
            m_sWGS = (D3DPOINT *) malloc(sizeof(D3DPOINT) * m_NumOfPts);
        } else {
            return nullptr;
        }
        
        
#pragma omp parallel for schedule(guided)
        for (int i = 0; i < m_NumOfPts; i++) {
            
            double x = _ps[i].m_X;
            double y = _ps[i].m_Y;
            
            double X = x - 500000;
            double Y = y;
            if (hemis < 0)
                Y = Y - 10000000;
            
            int S = ( ( _param.utm_zone* 6 ) - 183 );
            double lat =  Y / ( 6366197.724 * 0.9996 );
            double v = ( c / ( sqrt( 1.0 + ( e2cuadrada * ( cos(lat)*cos(lat) ) ) ) ) ) * 0.9996;
            double a = X / v;
            
            
            double a1 = sin( 2* lat );
            double a2 = a1* ( cos(lat)*cos(lat) );
            double j2 = lat + ( a1 / 2.0 );
            double j4 = ( ( 3 * j2 ) + a2 ) / 4.0;
            double j6 = ( ( 5 * j4 ) + ( a2 * ( cos(lat)*cos(lat) )) ) / 3;
            double alfa = ( 3 / 4.0 ) * e2cuadrada;
            double beta = ( 5 / 3.0 ) * alfa * alfa;
            double gama = ( 35 / 27.0 ) * alfa * alfa * alfa;
            double Bm = 0.9996 * c * ( lat - alfa * j2 + beta * j4 - gama * j6 );
            double b = ( Y - Bm ) / v;
            double Epsi = ( ( e2cuadrada * a *a ) / 2.0 ) * ( cos(lat)*cos(lat) );
            double Eps = a * ( 1 - ( Epsi / 3.0 ) );
            double nab = ( b * ( 1 - Epsi ) ) + lat;
            double senoheps = ( exp(Eps) - exp(-Eps) ) / 2;
            double Delt = atan(senoheps / (cos(nab) ) );
            double TaO = atan(cos(Delt) * tan(nab));
            double longitude = (Delt * (180/PI) ) + S;
            
            double latitude = ( lat + ( 1 + e2cuadrada * (cos(lat)*cos(lat)) - ( 3/2.0 )* e2cuadrada * sin(lat) * cos(lat) * ( TaO - lat ) )* ( TaO - lat ) ) * (180/PI);
            
            m_sWGS[i].m_Y = latitude;
            m_sWGS[i].m_X = longitude;
            m_sWGS[i].m_Z = _ps[i].m_Z;
        }
        return m_sWGS;
    }
    
}

D3DPOINT *ps2wgs_3D_vector(TransParam _param, int _numofpts, vector<D3DPOINT> &_ps)
{
    int m_NumOfPts = _numofpts;
    if(_param.projection == 1)
    {
        bool m_bHemisphere = _param.bHemisphere;
        double a = _param.a;
        double e = _param.e;
        double phi_c = _param.phi_c;
        double lambda_0 = _param.lambda_0;
        
        int pm = _param.pm;
        double t_c = _param.t_c;
        double m_c = _param.m_c;
            
        D3DPOINT *m_sWGS;
        
        if (m_NumOfPts > 0) {
            m_sWGS = (D3DPOINT *) malloc(sizeof(D3DPOINT) * m_NumOfPts);
        } else {
            return nullptr;
        }
        
        double e2 = e * e;
        double e4 = e2 * e2;
        double e6 = e2 * e4;
        double e8 = e4 * e4;
        
#pragma omp parallel for schedule(guided)
        for (int i = 0; i < m_NumOfPts; i++) {
          double x = _ps[i].m_X * pm;
          double y = _ps[i].m_Y * pm;
            
            double rho = sqrt(pow(x, 2) + pow(y, 2));
            double t = rho * t_c / (a * m_c);
            
            double chi = PI / 2 - 2 * atan(t);
            double phi = chi + (e2 / 2 + 5 * e4 / 24 + e6 / 12 + 13 * e8 / 360) * sin(2 * chi) + (7 * e4 / 48 + 29 * e6 / 240 + 811 * e8 / 11520) * sin(4 * chi) +
                (7 * e6 / 120 + 81 * e8 / 1120) * sin(6 * chi) + (4279 * e8 / 161280) * sin(8 * chi);
            
            double lambda = lambda_0 + atan2(x, -y);
            phi = pm * phi;
            lambda = pm * lambda;
            if (lambda > PI) {
                lambda = lambda - 2 * PI;
            } else if (lambda < -PI) {
                lambda = lambda + 2 * PI;
            }
            
            m_sWGS[i].m_Y = RadToDeg * phi;
            m_sWGS[i].m_X = RadToDeg * lambda;
            m_sWGS[i].m_Z = _ps[i].m_Z;
        }
            
        return m_sWGS;
    }
    else
    {
        int hemis = _param.pm; //1 = north, -1 = south
        double sa = _param.sa;
        double sb = _param.sb;
        double e2 = _param.e2;
        double e2cuadrada = _param.e2cuadrada;
        double c = _param.c;
            
        D3DPOINT *m_sWGS;
        
        if (m_NumOfPts > 0) {
            m_sWGS = (D3DPOINT *) malloc(sizeof(D3DPOINT) * m_NumOfPts);
        } else {
            return nullptr;
        }
        
        
#pragma omp parallel for schedule(guided)
        for (int i = 0; i < m_NumOfPts; i++) {
            
            double x = _ps[i].m_X;
            double y = _ps[i].m_Y;
            
            double X = x - 500000;
            double Y = y;
            if (hemis < 0)
                Y = Y - 10000000;
            
            int S = ( ( _param.utm_zone* 6 ) - 183 );
            double lat =  Y / ( 6366197.724 * 0.9996 );
            double v = ( c / ( sqrt( 1.0 + ( e2cuadrada * ( cos(lat)*cos(lat) ) ) ) ) ) * 0.9996;
            double a = X / v;
            
            
            double a1 = sin( 2* lat );
            double a2 = a1* ( cos(lat)*cos(lat) );
            double j2 = lat + ( a1 / 2.0 );
            double j4 = ( ( 3 * j2 ) + a2 ) / 4.0;
            double j6 = ( ( 5 * j4 ) + ( a2 * ( cos(lat)*cos(lat) )) ) / 3;
            double alfa = ( 3 / 4.0 ) * e2cuadrada;
            double beta = ( 5 / 3.0 ) * alfa * alfa;
            double gama = ( 35 / 27.0 ) * alfa * alfa * alfa;
            double Bm = 0.9996 * c * ( lat - alfa * j2 + beta * j4 - gama * j6 );
            double b = ( Y - Bm ) / v;
            double Epsi = ( ( e2cuadrada * a *a ) / 2.0 ) * ( cos(lat)*cos(lat) );
            double Eps = a * ( 1 - ( Epsi / 3.0 ) );
            double nab = ( b * ( 1 - Epsi ) ) + lat;
            double senoheps = ( exp(Eps) - exp(-Eps) ) / 2;
            double Delt = atan(senoheps / (cos(nab) ) );
            double TaO = atan(cos(Delt) * tan(nab));
            double longitude = (Delt * (180/PI) ) + S;
            
            double latitude = ( lat + ( 1 + e2cuadrada * (cos(lat)*cos(lat)) - ( 3/2.0 )* e2cuadrada * sin(lat) * cos(lat) * ( TaO - lat ) )* ( TaO - lat ) ) * (180/PI);
            
            m_sWGS[i].m_Y = latitude;
            m_sWGS[i].m_X = longitude;
            m_sWGS[i].m_Z = _ps[i].m_Z;
        }
        return m_sWGS;
    }
    
}

double GetVCPsIPsfromFRPCc(const double * const *rpc, const uint8 numofparam, const double *imageparam, CSize imagesize, vector<D3DPOINT> &VCPs, vector<D2DPOINT> &IPs)
{
    double mid_H;
    
    double VCP_interval = 10;
    double VCP_intervla_Z = 6;
    
    double minLon = (double) (-1.0 * rpc[1][2] + rpc[0][2]);
    double maxLon = (double) (1.0 * rpc[1][2] + rpc[0][2]);
    double minLat = (double) (-1.0 * rpc[1][3] + rpc[0][3]);
    double maxLat = (double) (1.0 * rpc[1][3] + rpc[0][3]);
    double minH = (double) (-1.0 * rpc[1][4] + rpc[0][4]);
    double maxH = (double) (1.0 * rpc[1][4] + rpc[0][4]);
    
    if(minH < 0)
        minH = 0;
    
    mid_H = minH + (maxH - minH)/2.0;
    
    double temp;
    if(minLon > maxLon)
    {
        temp = minLon;
        minLon = maxLon;
        maxLon = temp;
    }
    
    if(minLat > maxLat)
    {
        temp = minLat;
        minLat = maxLat;
        maxLat = temp;
    }
    
    //printf("minmax Lon %f\t%f\t%f\n",minLon,maxLon,rpc[0][2]);
    //printf("minmax Lat %f\t%f\t%f\n",minLat,maxLat,rpc[0][3]);
    //printf("minmax H %f\t%f\t%f\n",minH,maxH,rpc[0][4]);
    
    double X_dist = maxLon - minLon;
    double Y_dist = maxLat - minLat;
    double dist = (X_dist > Y_dist) ? Y_dist : X_dist;
    double Z_dist = maxH - minH;
    double X_interval = X_dist/VCP_interval;
    double Y_interval = Y_dist/VCP_interval;
    double Z_interval = Z_dist/VCP_intervla_Z;
    
    //printf("inverval %f\t%f\t%f\n",X_interval,Y_interval,Z_interval);
    for(int i=0;i<VCP_intervla_Z;i++) //Z
    {
        for(int j=0;j<VCP_interval;j++) //Y
        {
            for(int k=0;k<VCP_interval;k++) //X
            {
                D3DPOINT GP;
                GP.m_X = minLon + k*X_interval;
                GP.m_Y = minLat + j*Y_interval;
                GP.m_Z = minH + i*Z_interval;
                
                D2DPOINT IP = GetObjectToImageRPC_single_mpp(rpc,numofparam,imageparam,GP);
                
                //if(IP.m_X > 0 && IP.m_X < imagesize.width && IP.m_Y > 0 && IP.m_Y < imagesize.height)
                {
                    VCPs.push_back(GP);
                    IPs.push_back(IP);
                }
            }
        }
    }
    
    return mid_H;
}

double GetVCPsIPsfromFRPCc_Collinear(const double * const *rpc, TransParam param, CSize imagesize, EO eo, CAMERA_INFO camera, vector<D3DPOINT> &VCPs, vector<D2DPOINT> &IPs)
{
    double mid_H;
    
    double VCP_interval = 10;
    double VCP_intervla_Z = 6;
    
    double minLon = (double) (-1.0 * rpc[1][2] + rpc[0][2]);
    double maxLon = (double) (1.0 * rpc[1][2] + rpc[0][2]);
    double minLat = (double) (-1.0 * rpc[1][3] + rpc[0][3]);
    double maxLat = (double) (1.0 * rpc[1][3] + rpc[0][3]);
    double minH = (double) (-1.0 * rpc[1][4] + rpc[0][4]);
    double maxH = (double) (1.0 * rpc[1][4] + rpc[0][4]);
    
    if(minH < 0)
        minH = 0;
    mid_H = minH + (maxH - minH)/2.0;
    
    double temp;
    if(minLon > maxLon)
    {
        temp = minLon;
        minLon = maxLon;
        maxLon = temp;
    }
    
    if(minLat > maxLat)
    {
        temp = minLat;
        minLat = maxLat;
        maxLat = temp;
    }
    
    //printf("minmax Lon %f\t%f\t%f\n",minLon,maxLon,rpc[0][2]);
    //printf("minmax Lat %f\t%f\t%f\n",minLat,maxLat,rpc[0][3]);
    //printf("minmax H %f\t%f\t%f\n",minH,maxH,rpc[0][4]);
    
    double X_dist = maxLon - minLon;
    double Y_dist = maxLat - minLat;
    double dist = (X_dist > Y_dist) ? Y_dist : X_dist;
    double Z_dist = maxH - minH;
    double X_interval = X_dist/VCP_interval;
    double Y_interval = Y_dist/VCP_interval;
    double Z_interval = Z_dist/VCP_intervla_Z;
    
    //printf("inverval %f\t%f\t%f\n",X_interval,Y_interval,Z_interval);
    for(int i=0;i<VCP_intervla_Z;i++) //Z
    {
        for(int j=0;j<VCP_interval;j++) //Y
        {
            for(int k=0;k<VCP_interval;k++) //X
            {
                D3DPOINT GP;
                GP.m_X = minLon + k*X_interval;
                GP.m_Y = minLat + j*Y_interval;
                GP.m_Z = minH + i*Z_interval;
                
                D3DPOINT temppt = wgs2ps_single_3D(param,GP);
                
                D2DPOINT IP_photo = GetPhotoCoordinate_single(temppt, eo, camera, eo.m_Rm);
                D2DPOINT IP = PhotoToImage_single(IP_photo, camera.m_CCDSize, imagesize);
                //if(IP.m_X > 0 && IP.m_X < imagesize.width && IP.m_Y > 0 && IP.m_Y < imagesize.height)
                {
                    VCPs.push_back(GP);
                    IPs.push_back(IP);
                }
            }
        }
    }
    
    return mid_H;
}

double ** GetRPCsfromVCPsIPs(const double * const *rpc, const uint8 numofparam, const double *imageparam, vector<D3DPOINT> &VCPs, vector<D2DPOINT> &IPs)
{
    double ** IRPCs = (double**)calloc(7, sizeof(double*));
    IRPCs[0] = (double*)calloc(5, sizeof(double)); //lineoffset, sampleoffset, longoffset, latoffset, heightoffset
    IRPCs[1] = (double*)calloc(5, sizeof(double)); //linescale, samplescale, longscale, latscale, heightscale
    IRPCs[2] = (double*)calloc(20, sizeof(double)); //linenumcoef
    IRPCs[3] = (double*)calloc(20, sizeof(double)); //linedencoef
    IRPCs[4] = (double*)calloc(20, sizeof(double)); //samplenumcoef
    IRPCs[5] = (double*)calloc(20, sizeof(double)); //sampledencoef
    IRPCs[6] = (double*)calloc(2, sizeof(double)); //errbias, errand for worldview
    
    for(int i=0;i<5;i++)
    {
        IRPCs[0][i] = rpc[0][i];
        IRPCs[1][i] = rpc[1][i];
    }

    for(int i=0;i<2;i++)
        IRPCs[6][i] = rpc[6][i];
    
    IRPCs[3][0] = 1.0;
    IRPCs[5][0] = 1.0;
    
    int count_GCPs = VCPs.size();
    
    for(int whole_iter = 0 ; whole_iter < 2 ; whole_iter++)
    {
        GMA_double *M_matrix = GMA_double_create(count_GCPs,39);
        GMA_double *MT_matrix = GMA_double_create(39,count_GCPs);
        GMA_double *MTM_matrix = GMA_double_create(39,39);
        GMA_double *MTMI_matrix = GMA_double_create(39,39);
        
        GMA_double *R_matrix = GMA_double_create(count_GCPs,1);
        GMA_double *MTR_matrix = GMA_double_create(39,1);
        GMA_double *J_matrix = GMA_double_create(39,1);
        
        //printf("done matrix create 1\n");
        
        //L => X(c), P => Y(r), H => Z
        double L,P,H,r,c;
        for(int i=0;i<count_GCPs;i++)
        {
            L       = (VCPs[i].m_X - rpc[0][2])/rpc[1][2];
            P       = (VCPs[i].m_Y - rpc[0][3])/rpc[1][3];
            H       = (VCPs[i].m_Z - rpc[0][4])/rpc[1][4];
            
            r       = (IPs[i].m_Y - imageparam[0] - rpc[0][0])/rpc[1][0]; //Line
            c       = (IPs[i].m_X - imageparam[1] - rpc[0][1])/rpc[1][1]; //sample
            
            
            double temp = r;
            if(whole_iter == 1)
                temp = c;
            
            r       = P;
            c       = L;
            P       = temp;
            
            /*
            if(whole_iter == 1)
                P = temp;
            */
            M_matrix->val[i][0] = 1.0;
            M_matrix->val[i][1] = c;
            M_matrix->val[i][2] = r;
            
            M_matrix->val[i][3] = H;
            M_matrix->val[i][4] = c*r;
            M_matrix->val[i][5] = c*H;
            
            M_matrix->val[i][6] = r*H;
            M_matrix->val[i][7] = c*c;
            M_matrix->val[i][8] = r*r;
            
            M_matrix->val[i][9] = H*H;
            M_matrix->val[i][10] = c*r*H;
            M_matrix->val[i][11] = c*c*c;
            
            M_matrix->val[i][12] = c*r*r;
            M_matrix->val[i][13] = c*H*H;
            M_matrix->val[i][14] = c*c*r;
            
            M_matrix->val[i][15] = r*r*r;
            M_matrix->val[i][16] = r*H*H;
            M_matrix->val[i][17] = c*c*H;
            
            M_matrix->val[i][18] = r*r*H;
            M_matrix->val[i][19] = H*H*H;
            
            
            M_matrix->val[i][20] = -P*c;
            M_matrix->val[i][21] = -P*r;
            
            M_matrix->val[i][22] = -P*H;
            M_matrix->val[i][23] = -P*c*r;
            M_matrix->val[i][24] = -P*c*H;
            
            M_matrix->val[i][25] = -P*r*H;
            M_matrix->val[i][26] = -P*c*c;
            M_matrix->val[i][27] = -P*r*r;
            
            M_matrix->val[i][28] = -P*H*H;
            M_matrix->val[i][29] = -P*c*r*H;
            M_matrix->val[i][30] = -P*c*c*c;
            
            M_matrix->val[i][31] = -P*c*r*r;
            M_matrix->val[i][32] = -P*c*H*H;
            M_matrix->val[i][33] = -P*c*c*r;
            
            M_matrix->val[i][34] = -P*r*r*r;
            M_matrix->val[i][35] = -P*r*H*H;
            M_matrix->val[i][36] = -P*c*c*H;
            
            M_matrix->val[i][37] = -P*r*r*H;
            M_matrix->val[i][38] = -P*H*H*H;
            
            
            R_matrix->val[i][0] = P;
        }
        
        //GMA_double_printf(M_matrix);
        //GMA_double_printf(R_matrix);
        
        //printf("done matrix M R assign\n");
        
        GMA_double_Tran(M_matrix,MT_matrix);
        GMA_double_mul(MT_matrix,M_matrix,MTM_matrix);
        GMA_double_inv(MTM_matrix,MTMI_matrix);
        GMA_double_mul(MT_matrix,R_matrix,MTR_matrix);
        GMA_double_mul(MTMI_matrix,MTR_matrix,J_matrix);
        
        //printf("done matrix cal 1\n");
        
        GMA_double_destroy(MTM_matrix);
        GMA_double_destroy(MTMI_matrix);
        GMA_double_destroy(MTR_matrix);
        
        double sigma_pre = 100000;
        double sigma_ratio = 100000;
        int max_iter = 10;
        int iteration = 0;
        
        while(iteration < max_iter && sigma_ratio > 0.1)
        {
            GMA_double *DEN_matrix = GMA_double_create(20,1);
            GMA_double *W_matrix = GMA_double_create(count_GCPs,count_GCPs);
            
            DEN_matrix->val[0][0] = 1.0;
            for(int i=1;i<20;i++)
            {
                DEN_matrix->val[i][0] = J_matrix->val[i+19][0];
                //printf("coef %d\t %Lf\n",i,DEN_matrix->val[i][0] );
            }
            
            //GMA_double_printf(DEN_matrix);
            //GMA_double_printf(J_matrix);
            
            //GMA_double_destroy(J_matrix);
            
            //printf("done matrix create 2\n");
            
            //weight matrix
            for(int i=0;i<count_GCPs;i++)
            {
                GMA_double *B_matrix = GMA_double_create(1,20);
                GMA_double *BDEN_matrix = GMA_double_create(1,1);
          
                L       = (VCPs[i].m_X - rpc[0][2])/rpc[1][2];
                P       = (VCPs[i].m_Y - rpc[0][3])/rpc[1][3];
                H       = (VCPs[i].m_Z - rpc[0][4])/rpc[1][4];
                
                r       = (IPs[i].m_Y - imageparam[0] - rpc[0][0])/rpc[1][0]; //Line
                c       = (IPs[i].m_X - imageparam[1] - rpc[0][1])/rpc[1][1]; //sample
               
                
                r       = P;
                c       = L;
                
                
                B_matrix->val[0][0] = 1.0;
                B_matrix->val[0][1] = c;
                B_matrix->val[0][2] = r;
                
                B_matrix->val[0][3] = H;
                B_matrix->val[0][4] = c*r;
                B_matrix->val[0][5] = c*H;
                
                B_matrix->val[0][6] = r*H;
                B_matrix->val[0][7] = c*c;
                B_matrix->val[0][8] = r*r;
                
                B_matrix->val[0][9] = H*H;
                B_matrix->val[0][10] = c*r*H;
                B_matrix->val[0][11] = c*c*c;
                
                B_matrix->val[0][12] = c*r*r;
                B_matrix->val[0][13] = c*H*H;
                B_matrix->val[0][14] = c*c*r;
                
                B_matrix->val[0][15] = r*r*r;
                B_matrix->val[0][16] = r*H*H;
                B_matrix->val[0][17] = c*c*H;
                
                B_matrix->val[0][18] = r*r*H;
                B_matrix->val[0][19] = H*H*H;
                
                GMA_double_mul(B_matrix,DEN_matrix,BDEN_matrix);
                double B = BDEN_matrix->val[0][0];
                
                W_matrix->val[i][i] = (double)(1.0/(B*B));
                
                GMA_double_destroy(B_matrix);
                GMA_double_destroy(BDEN_matrix);
            }
            
            //printf("done matrix W assign\n");
            //GMA_double_printf(W_matrix);
            
            GMA_double_destroy(DEN_matrix);
            
            GMA_double *MTW_matrix = GMA_double_create(39,count_GCPs);
            GMA_double *MTWM_matrix = GMA_double_create(39,39);
            GMA_double *MTWMI_matrix = GMA_double_create(39,39);
            GMA_double *MTWR_matrix = GMA_double_create(39,1);
            //GMA_double *JJ_matrix = GMA_double_create(39,1);
            
            GMA_double *WM_matrix = GMA_double_create(count_GCPs,39);
            GMA_double *WMJ_matrix = GMA_double_create(count_GCPs,1);
            GMA_double *WR_matrix = GMA_double_create(count_GCPs,1);
            GMA_double *V_matrix = GMA_double_create(count_GCPs,1);
            
            //printf("done matrix create 3\n");
            
            GMA_double_mul(MT_matrix,W_matrix,MTW_matrix);
            //printf("MTW_matrix\n");
            GMA_double_mul(MTW_matrix,M_matrix,MTWM_matrix);
            //printf("MTWM_matrix\n");
            GMA_double_inv(MTWM_matrix,MTWMI_matrix);
            //printf("MTWMI_matrix\n");
            GMA_double_mul(MTW_matrix,R_matrix,MTWR_matrix);
            //printf("MTWR_matrix\n");
            GMA_double_mul(MTWMI_matrix,MTWR_matrix,J_matrix);
            
            //printf("done matrix cal 2 1\n");
            
            GMA_double_mul(W_matrix,M_matrix,WM_matrix);
            GMA_double_mul(W_matrix,R_matrix,WR_matrix);
            GMA_double_sub(WMJ_matrix,WR_matrix,V_matrix);
            
            //printf("done matrix cal 2 2\n");
            
            
            
            GMA_double_destroy(W_matrix);
            
            GMA_double_destroy(MTW_matrix);
            GMA_double_destroy(MTWM_matrix);
            GMA_double_destroy(MTWMI_matrix);
            GMA_double_destroy(MTWR_matrix);
            
            GMA_double_destroy(WM_matrix);
            GMA_double_destroy(WMJ_matrix);
            GMA_double_destroy(WR_matrix);
            
            double sum=0;
            for(int i=0;i<count_GCPs;i++)
                sum += (V_matrix->val[i][0])*(V_matrix->val[i][0]);
            
            GMA_double_destroy(V_matrix);
            
            double sigma = sqrt(sum/(count_GCPs-39));
            
            sigma_ratio = fabs(sigma_pre - sigma)/sigma_pre;
            sigma_pre = sigma;
            
            //printf("iteration %d\tsigma %f\t%f\n", iteration, sigma,sigma_ratio);
            
            //GMA_double_printf(J_matrix);
            
            iteration++;
            
        }
        
        if(whole_iter == 1)
        {
            for(int i=0;i<20;i++)
            {
                IRPCs[4][i] = J_matrix->val[i][0];
                //printf("rpc %d\t%e\n",i+1,IRPCs[4][i]);
            }
            
            IRPCs[5][0] = 1.0;
            //printf("rpc 21\t%e\n",IRPCs[5][0]);
            for(int i=1;i<20;i++)
            {
                IRPCs[5][i] = J_matrix->val[i+19][0];
                //printf("rpc %d\t%e\n",i+21,IRPCs[5][i]);
            }
        }
        else
        {
            for(int i=0;i<20;i++)
            {
                IRPCs[2][i] = J_matrix->val[i][0];
                //printf("rpc %d\t%e\n",i+1,IRPCs[2][i]);
            }
            
            IRPCs[3][0] = 1.0;
            //printf("rpc 21\t%e\n",IRPCs[3][0]);
            for(int i=1;i<20;i++)
            {
                IRPCs[3][i] = J_matrix->val[i+19][0];
                //printf("rpc %d\t%e\n",i+21,IRPCs[3][i]);
            }
        }
        
        GMA_double_destroy(M_matrix);
        GMA_double_destroy(MT_matrix);
        GMA_double_destroy(R_matrix);
        
        GMA_double_destroy(J_matrix);
    }
    
    return IRPCs;
}

double ** GetIRPCsfromVCPsIPs(const double * const *rpc, const uint8 numofparam, const double *imageparam, vector<D3DPOINT> &VCPs, vector<D2DPOINT> &IPs)
{
    double ** IRPCs = (double**)calloc(7, sizeof(double*));
    IRPCs[0] = (double*)calloc(5, sizeof(double)); //lineoffset, sampleoffset, longoffset, latoffset, heightoffset
    IRPCs[1] = (double*)calloc(5, sizeof(double)); //linescale, samplescale, longscale, latscale, heightscale
    IRPCs[2] = (double*)calloc(20, sizeof(double)); //linenumcoef
    IRPCs[3] = (double*)calloc(20, sizeof(double)); //linedencoef
    IRPCs[4] = (double*)calloc(20, sizeof(double)); //samplenumcoef
    IRPCs[5] = (double*)calloc(20, sizeof(double)); //sampledencoef
    IRPCs[6] = (double*)calloc(2, sizeof(double)); //errbias, errand for worldview
    
    for(int i=0;i<5;i++)
    {
        IRPCs[0][i] = rpc[0][i];
        IRPCs[1][i] = rpc[1][i];
    }

    for(int i=0;i<2;i++)
        IRPCs[6][i] = rpc[6][i];
    
    IRPCs[3][0] = 1.0;
    IRPCs[5][0] = 1.0;
    
    int count_GCPs = VCPs.size();
    
    for(int whole_iter = 0 ; whole_iter < 2 ; whole_iter++)
    {
        GMA_double *M_matrix = GMA_double_create(count_GCPs,39);
        GMA_double *MT_matrix = GMA_double_create(39,count_GCPs);
        GMA_double *MTM_matrix = GMA_double_create(39,39);
        GMA_double *MTMI_matrix = GMA_double_create(39,39);
        
        GMA_double *R_matrix = GMA_double_create(count_GCPs,1);
        GMA_double *MTR_matrix = GMA_double_create(39,1);
        GMA_double *J_matrix = GMA_double_create(39,1);
        
        //printf("done matrix create 1\n");
        
        //L => X(c), P => Y(r), H => Z
        double L,P,H,r,c;
        for(int i=0;i<count_GCPs;i++)
        {
            L       = (VCPs[i].m_X - rpc[0][2])/rpc[1][2];
            P       = (VCPs[i].m_Y - rpc[0][3])/rpc[1][3];
            H       = (VCPs[i].m_Z - rpc[0][4])/rpc[1][4];
            
            r       = (IPs[i].m_Y - imageparam[0] - rpc[0][0])/rpc[1][0]; //Line
            c       = (IPs[i].m_X - imageparam[1] - rpc[0][1])/rpc[1][1]; //sample
            
            /*
            double temp = r;
            if(whole_iter == 1)
                temp = c;
            
            r       = P;
            c       = L;
            P       = temp;
            */
           
            if(whole_iter == 1)
                P = L;
            
            M_matrix->val[i][0] = 1.0;
            M_matrix->val[i][1] = c;
            M_matrix->val[i][2] = r;
            
            M_matrix->val[i][3] = H;
            M_matrix->val[i][4] = c*r;
            M_matrix->val[i][5] = c*H;
            
            M_matrix->val[i][6] = r*H;
            M_matrix->val[i][7] = c*c;
            M_matrix->val[i][8] = r*r;
            
            M_matrix->val[i][9] = H*H;
            M_matrix->val[i][10] = c*r*H;
            M_matrix->val[i][11] = c*c*c;
            
            M_matrix->val[i][12] = c*r*r;
            M_matrix->val[i][13] = c*H*H;
            M_matrix->val[i][14] = c*c*r;
            
            M_matrix->val[i][15] = r*r*r;
            M_matrix->val[i][16] = r*H*H;
            M_matrix->val[i][17] = c*c*H;
            
            M_matrix->val[i][18] = r*r*H;
            M_matrix->val[i][19] = H*H*H;
            
            
            M_matrix->val[i][20] = -P*c;
            M_matrix->val[i][21] = -P*r;
            
            M_matrix->val[i][22] = -P*H;
            M_matrix->val[i][23] = -P*c*r;
            M_matrix->val[i][24] = -P*c*H;
            
            M_matrix->val[i][25] = -P*r*H;
            M_matrix->val[i][26] = -P*c*c;
            M_matrix->val[i][27] = -P*r*r;
            
            M_matrix->val[i][28] = -P*H*H;
            M_matrix->val[i][29] = -P*c*r*H;
            M_matrix->val[i][30] = -P*c*c*c;
            
            M_matrix->val[i][31] = -P*c*r*r;
            M_matrix->val[i][32] = -P*c*H*H;
            M_matrix->val[i][33] = -P*c*c*r;
            
            M_matrix->val[i][34] = -P*r*r*r;
            M_matrix->val[i][35] = -P*r*H*H;
            M_matrix->val[i][36] = -P*c*c*H;
            
            M_matrix->val[i][37] = -P*r*r*H;
            M_matrix->val[i][38] = -P*H*H*H;
            
            
            R_matrix->val[i][0] = P;
        }
        
        //GMA_double_printf(M_matrix);
        //GMA_double_printf(R_matrix);
        
        //printf("done matrix M R assign\n");
        
        GMA_double_Tran(M_matrix,MT_matrix);
        GMA_double_mul(MT_matrix,M_matrix,MTM_matrix);
        GMA_double_inv(MTM_matrix,MTMI_matrix);
        GMA_double_mul(MT_matrix,R_matrix,MTR_matrix);
        GMA_double_mul(MTMI_matrix,MTR_matrix,J_matrix);
        
        //printf("done matrix cal 1\n");
        
        GMA_double_destroy(MTM_matrix);
        GMA_double_destroy(MTMI_matrix);
        GMA_double_destroy(MTR_matrix);
        
        double sigma_pre = 100000;
        double sigma_ratio = 100000;
        int max_iter = 10;
        int iteration = 0;
        
        while(iteration < max_iter && sigma_ratio > 0.1)
        {
            GMA_double *DEN_matrix = GMA_double_create(20,1);
            GMA_double *W_matrix = GMA_double_create(count_GCPs,count_GCPs);
            
            DEN_matrix->val[0][0] = 1.0;
            for(int i=1;i<20;i++)
            {
                DEN_matrix->val[i][0] = J_matrix->val[i+19][0];
                //printf("coef %d\t %Lf\n",i,DEN_matrix->val[i][0] );
            }
            
            //GMA_double_printf(DEN_matrix);
            //GMA_double_printf(J_matrix);
            
            //GMA_double_destroy(J_matrix);
            
            //printf("done matrix create 2\n");
            
            //weight matrix
            for(int i=0;i<count_GCPs;i++)
            {
                GMA_double *B_matrix = GMA_double_create(1,20);
                GMA_double *BDEN_matrix = GMA_double_create(1,1);
          
                L       = (VCPs[i].m_X - rpc[0][2])/rpc[1][2];
                P       = (VCPs[i].m_Y - rpc[0][3])/rpc[1][3];
                H       = (VCPs[i].m_Z - rpc[0][4])/rpc[1][4];
                
                r       = (IPs[i].m_Y - imageparam[0] - rpc[0][0])/rpc[1][0]; //Line
                c       = (IPs[i].m_X - imageparam[1] - rpc[0][1])/rpc[1][1]; //sample
               
                /*
                r       = P;
                c       = L;
                */
                
                B_matrix->val[0][0] = 1.0;
                B_matrix->val[0][1] = c;
                B_matrix->val[0][2] = r;
                
                B_matrix->val[0][3] = H;
                B_matrix->val[0][4] = c*r;
                B_matrix->val[0][5] = c*H;
                
                B_matrix->val[0][6] = r*H;
                B_matrix->val[0][7] = c*c;
                B_matrix->val[0][8] = r*r;
                
                B_matrix->val[0][9] = H*H;
                B_matrix->val[0][10] = c*r*H;
                B_matrix->val[0][11] = c*c*c;
                
                B_matrix->val[0][12] = c*r*r;
                B_matrix->val[0][13] = c*H*H;
                B_matrix->val[0][14] = c*c*r;
                
                B_matrix->val[0][15] = r*r*r;
                B_matrix->val[0][16] = r*H*H;
                B_matrix->val[0][17] = c*c*H;
                
                B_matrix->val[0][18] = r*r*H;
                B_matrix->val[0][19] = H*H*H;
                
                GMA_double_mul(B_matrix,DEN_matrix,BDEN_matrix);
                double B = BDEN_matrix->val[0][0];
                
                W_matrix->val[i][i] = (double)(1.0/(B*B));
                
                GMA_double_destroy(B_matrix);
                GMA_double_destroy(BDEN_matrix);
            }
            
            //printf("done matrix W assign\n");
            //GMA_double_printf(W_matrix);
            
            GMA_double_destroy(DEN_matrix);
            
            GMA_double *MTW_matrix = GMA_double_create(39,count_GCPs);
            GMA_double *MTWM_matrix = GMA_double_create(39,39);
            GMA_double *MTWMI_matrix = GMA_double_create(39,39);
            GMA_double *MTWR_matrix = GMA_double_create(39,1);
            //GMA_double *JJ_matrix = GMA_double_create(39,1);
            
            GMA_double *WM_matrix = GMA_double_create(count_GCPs,39);
            GMA_double *WMJ_matrix = GMA_double_create(count_GCPs,1);
            GMA_double *WR_matrix = GMA_double_create(count_GCPs,1);
            GMA_double *V_matrix = GMA_double_create(count_GCPs,1);
            
            //printf("done matrix create 3\n");
            
            GMA_double_mul(MT_matrix,W_matrix,MTW_matrix);
            //printf("MTW_matrix\n");
            GMA_double_mul(MTW_matrix,M_matrix,MTWM_matrix);
            //printf("MTWM_matrix\n");
            GMA_double_inv(MTWM_matrix,MTWMI_matrix);
            //printf("MTWMI_matrix\n");
            GMA_double_mul(MTW_matrix,R_matrix,MTWR_matrix);
            //printf("MTWR_matrix\n");
            GMA_double_mul(MTWMI_matrix,MTWR_matrix,J_matrix);
            
            //printf("done matrix cal 2 1\n");
            
            GMA_double_mul(W_matrix,M_matrix,WM_matrix);
            GMA_double_mul(W_matrix,R_matrix,WR_matrix);
            GMA_double_sub(WMJ_matrix,WR_matrix,V_matrix);
            
            //printf("done matrix cal 2 2\n");
            
            
            
            GMA_double_destroy(W_matrix);
            
            GMA_double_destroy(MTW_matrix);
            GMA_double_destroy(MTWM_matrix);
            GMA_double_destroy(MTWMI_matrix);
            GMA_double_destroy(MTWR_matrix);
            
            GMA_double_destroy(WM_matrix);
            GMA_double_destroy(WMJ_matrix);
            GMA_double_destroy(WR_matrix);
            
            double sum=0;
            for(int i=0;i<count_GCPs;i++)
                sum += (V_matrix->val[i][0])*(V_matrix->val[i][0]);
            
            GMA_double_destroy(V_matrix);
            
            double sigma = sqrt(sum/(count_GCPs-39));
            
            sigma_ratio = fabs(sigma_pre - sigma)/sigma_pre;
            sigma_pre = sigma;
            
            //printf("iteration %d\tsigma %f\t%f\n", iteration, sigma,sigma_ratio);
            
            //GMA_double_printf(J_matrix);
            
            iteration++;
            
        }
        
        if(whole_iter == 1)
        {
            for(int i=0;i<20;i++)
            {
                IRPCs[4][i] = J_matrix->val[i][0];
                //printf("rpc %d\t%e\n",i+1,IRPCs[4][i]);
            }
            
            IRPCs[5][0] = 1.0;
            //printf("rpc 21\t%e\n",IRPCs[5][0]);
            for(int i=1;i<20;i++)
            {
                IRPCs[5][i] = J_matrix->val[i+19][0];
                //printf("rpc %d\t%e\n",i+21,IRPCs[5][i]);
            }
        }
        else
        {
            for(int i=0;i<20;i++)
            {
                IRPCs[2][i] = J_matrix->val[i][0];
                //printf("rpc %d\t%e\n",i+1,IRPCs[2][i]);
            }
            
            IRPCs[3][0] = 1.0;
            //printf("rpc 21\t%e\n",IRPCs[3][0]);
            for(int i=1;i<20;i++)
            {
                IRPCs[3][i] = J_matrix->val[i+19][0];
                //printf("rpc %d\t%e\n",i+21,IRPCs[3][i]);
            }
        }
        
        GMA_double_destroy(M_matrix);
        GMA_double_destroy(MT_matrix);
        GMA_double_destroy(R_matrix);
        
        GMA_double_destroy(J_matrix);
    }
    
    return IRPCs;
}

D2DPOINT* GetObjectToImageRPC(const double * const *_rpc, const uint8 _numofparam, const double *_imageparam, const uint16 _numofpts, D3DPOINT *_GP)
{
    D2DPOINT *IP;
    
    IP      = (D2DPOINT*)malloc(sizeof(D2DPOINT)*_numofpts);
    
#pragma omp parallel for schedule(guided)
    for(int i=0;i<_numofpts;i++)
    {
        double L, P, H, Line, Samp, deltaP, deltaR;
        double Coeff[4];

        L       = (_GP[i].m_X - _rpc[0][2])/_rpc[1][2];
        P       = (_GP[i].m_Y - _rpc[0][3])/_rpc[1][3];
        H       = (_GP[i].m_Z - _rpc[0][4])/_rpc[1][4];
        
        //printf("original %f\t%f\t%f\tL P H %f\t%f\t%f\n",_GP[i].m_X,_GP[i].m_Y,_GP[i].m_Z,L,P,H);
        
        if(L < -10.0 || L > 10.0)
        {
            if(_GP[i].m_X > 0)
                _GP[i].m_X = _GP[i].m_X - 360;
            else
                _GP[i].m_X = _GP[i].m_X + 360;
            
            L       = (_GP[i].m_X - _rpc[0][2])/_rpc[1][2];
            
            //printf("original %f\t%f\t%f\tL P H %f\t%f\t%f\n",_GP[i].m_X,_GP[i].m_Y,_GP[i].m_Z,L,P,H);
            
        }
        if(P < -10.0 || P > 10.0)
        {
            if(_GP[i].m_Y > 0)
                _GP[i].m_Y = _GP[i].m_Y - 360;
            else
                _GP[i].m_Y = _GP[i].m_Y + 360;
            
            P       = (_GP[i].m_Y - _rpc[0][3])/_rpc[1][3];
            
            //printf("original %f\t%f\t%f\tL P H %f\t%f\t%f\n",_GP[i].m_X,_GP[i].m_Y,_GP[i].m_Z,L,P,H);
        }
        
        //printf("L P H %f\t%f\t%f\n",L,P,H);
        
        for(int j=0;j<4;j++)
        {
            Coeff[j]    = _rpc[j+2][0]*1.0          + _rpc[j+2][1]*L            + _rpc[j+2][2]*P
                + _rpc[j+2][3]*H            + _rpc[j+2][4]*L*P          + _rpc[j+2][5]*L*H
                + _rpc[j+2][6]*P*H          + _rpc[j+2][7]*L*L          + _rpc[j+2][8]*P*P
                + _rpc[j+2][9]*H*H          + _rpc[j+2][10]*(P*L)*H     + _rpc[j+2][11]*(L*L)*L
                + _rpc[j+2][12]*(L*P)*P     + _rpc[j+2][13]*(L*H)*H     + _rpc[j+2][14]*(L*L)*P
                + _rpc[j+2][15]*(P*P)*P     + _rpc[j+2][16]*(P*H)*H     + _rpc[j+2][17]*(L*L)*H
                + _rpc[j+2][18]*(P*P)*H     + _rpc[j+2][19]*(H*H)*H;
        }

        Line     = ((Coeff[0]/Coeff[1])*_rpc[1][0] + _rpc[0][0]); //Line
        Samp     = ((Coeff[2]/Coeff[3])*_rpc[1][1] + _rpc[0][1]); //Sample
        deltaP      = _imageparam[0];
        deltaR      = _imageparam[1];

        IP[i].m_Y       = deltaP + Line;
        IP[i].m_X       = deltaR + Samp;
        
        if(IP[i].m_Y < 0)
            IP[i].m_Y = 0;
        
        if(IP[i].m_Y > _rpc[0][0] + _rpc[1][0]*1.2)
            IP[i].m_Y = _rpc[0][0] + _rpc[1][0]*1.2;
        
        if(IP[i].m_X < 0)
            IP[i].m_X = 0;
        
        if(IP[i].m_X > _rpc[0][1] + _rpc[1][1]*1.2)
            IP[i].m_X = _rpc[0][1] + _rpc[1][1]*1.2;
    }

    return IP;
}

D3DPOINT* GetImageHToObjectIRPC(const double * const *_rpc, const uint8 _numofparam, const double *_imageparam, vector<D3DPOINT> &_GP, vector<D2DPOINT> &_IPs)
{
    D3DPOINT *Object;
    
    long _numofpts = _GP.size();
    
    Object      = (D3DPOINT*)malloc(sizeof(D3DPOINT)*_numofpts);
    
#pragma omp parallel for schedule(guided)
    for(int i=0;i<_numofpts;i++)
    {
        double L, P, H, r, c;
        double Coeff[4];

        r       = (_IPs[i].m_Y - _imageparam[0] - _rpc[0][0])/_rpc[1][0]; //Line
        c       = (_IPs[i].m_X - _imageparam[1] - _rpc[0][1])/_rpc[1][1]; //sample
        H       = (_GP[i].m_Z - _rpc[0][4])/_rpc[1][4];
        
        P       = r;
        L       = c;
        //printf("original %f\t%f\t%f\tL P H %f\t%f\t%f\n",_GP[i].m_X,_GP[i].m_Y,_GP[i].m_Z,L,P,H);
        
        for(int j=0;j<4;j++)
        {
            Coeff[j]    = _rpc[j+2][0]*1.0          + _rpc[j+2][1]*L            + _rpc[j+2][2]*P
                + _rpc[j+2][3]*H            + _rpc[j+2][4]*L*P          + _rpc[j+2][5]*L*H
                + _rpc[j+2][6]*P*H          + _rpc[j+2][7]*L*L          + _rpc[j+2][8]*P*P
                + _rpc[j+2][9]*H*H          + _rpc[j+2][10]*(P*L)*H     + _rpc[j+2][11]*(L*L)*L
                + _rpc[j+2][12]*(L*P)*P     + _rpc[j+2][13]*(L*H)*H     + _rpc[j+2][14]*(L*L)*P
                + _rpc[j+2][15]*(P*P)*P     + _rpc[j+2][16]*(P*H)*H     + _rpc[j+2][17]*(L*L)*H
                + _rpc[j+2][18]*(P*P)*H     + _rpc[j+2][19]*(H*H)*H;
        }

        P       = Coeff[0]/Coeff[1];
        L       = Coeff[2]/Coeff[3];
        
        Object[i].m_Y = P*_rpc[1][3] + _rpc[0][3]; //Y
        Object[i].m_X = L*_rpc[1][2] + _rpc[0][2]; //X
        Object[i].m_Z = _GP[i].m_Z;//Z
        
        if(L < -10.0 || L > 10.0)
        {
            if(Object[i].m_X > 0)
                Object[i].m_X = Object[i].m_X - 360;
            else
                Object[i].m_X = Object[i].m_X + 360;
        }
        
        if(P < -10.0 || P > 10.0)
        {
            if(Object[i].m_Y > 0)
                Object[i].m_Y = Object[i].m_Y - 360;
            else
                Object[i].m_Y = Object[i].m_Y + 360;
            
        }
    }

    return Object;
}

D2DPOINT GetObjectToImageRPC_single2(const double * const *_rpc, const uint8 _numofparam, const double *_imageparam, D3DPOINT _GP)
{
    D2DPOINT IP;
    
    {
        double L, P, H, Line, Samp, deltaP, deltaR;
        double Coeff[4];

        L       = (_GP.m_X - _rpc[0][2])/_rpc[1][2];
        P       = (_GP.m_Y - _rpc[0][3])/_rpc[1][3];
        H       = (_GP.m_Z - _rpc[0][4])/_rpc[1][4];
        
        //printf("original %f\t%f\t%f\tL P H %f\t%f\t%f\n",_GP[i].m_X,_GP[i].m_Y,_GP[i].m_Z,L,P,H);
        
        if(L < -10.0 || L > 10.0)
        {
            if(_GP.m_X > 0)
                _GP.m_X = _GP.m_X - 360;
            else
                _GP.m_X = _GP.m_X + 360;
            
            L       = (_GP.m_X - _rpc[0][2])/_rpc[1][2];
            
            //printf("original %f\t%f\t%f\tL P H %f\t%f\t%f\n",_GP[i].m_X,_GP[i].m_Y,_GP[i].m_Z,L,P,H);
            
        }
        if(P < -10.0 || P > 10.0)
        {
            if(_GP.m_Y > 0)
                _GP.m_Y = _GP.m_Y - 360;
            else
                _GP.m_Y = _GP.m_Y + 360;
            
            P       = (_GP.m_Y - _rpc[0][3])/_rpc[1][3];
            
            //printf("original %f\t%f\t%f\tL P H %f\t%f\t%f\n",_GP[i].m_X,_GP[i].m_Y,_GP[i].m_Z,L,P,H);
        }
        
        //printf("L P H %f\t%f\t%f\n",L,P,H);
        
        for(int j=0;j<4;j++)
        {
            Coeff[j]    = _rpc[j+2][0]*1.0          + _rpc[j+2][1]*L            + _rpc[j+2][2]*P
                + _rpc[j+2][3]*H            + _rpc[j+2][4]*L*P          + _rpc[j+2][5]*L*H
                + _rpc[j+2][6]*P*H          + _rpc[j+2][7]*L*L          + _rpc[j+2][8]*P*P
                + _rpc[j+2][9]*H*H          + _rpc[j+2][10]*(P*L)*H     + _rpc[j+2][11]*(L*L)*L
                + _rpc[j+2][12]*(L*P)*P     + _rpc[j+2][13]*(L*H)*H     + _rpc[j+2][14]*(L*L)*P
                + _rpc[j+2][15]*(P*P)*P     + _rpc[j+2][16]*(P*H)*H     + _rpc[j+2][17]*(L*L)*H
                + _rpc[j+2][18]*(P*P)*H     + _rpc[j+2][19]*(H*H)*H;
        }

        Line     = ((Coeff[0]/Coeff[1])*_rpc[1][0] + _rpc[0][0]); //Line
        Samp     = ((Coeff[2]/Coeff[3])*_rpc[1][1] + _rpc[0][1]); //Sample
        deltaP      = _imageparam[0];
        deltaR      = _imageparam[1];

        IP.m_Y       = deltaP + Line;
        IP.m_X       = deltaR + Samp;
        
    }

    return IP;
}

D3DPOINT GetImageHToObjectIRPC_single(const double * const *_rpc, const uint8 _numofparam, const double *_imageparam, double GP_H, D2DPOINT _IPs)
{
    D3DPOINT Object;
    
    {
        double L, P, H, r, c;
        double Coeff[4];

        r       = (_IPs.m_Y - _imageparam[0] - _rpc[0][0])/_rpc[1][0]; //Line
        c       = (_IPs.m_X - _imageparam[1] - _rpc[0][1])/_rpc[1][1]; //sample
        H       = (GP_H - _rpc[0][4])/_rpc[1][4];
        
        P       = r;
        L       = c;
        //printf("original %f\t%f\t%f\tL P H %f\t%f\t%f\n",_GP[i].m_X,_GP[i].m_Y,_GP[i].m_Z,L,P,H);
        
        for(int j=0;j<4;j++)
        {
            Coeff[j]    = _rpc[j+2][0]*1.0          + _rpc[j+2][1]*L            + _rpc[j+2][2]*P
                + _rpc[j+2][3]*H            + _rpc[j+2][4]*L*P          + _rpc[j+2][5]*L*H
                + _rpc[j+2][6]*P*H          + _rpc[j+2][7]*L*L          + _rpc[j+2][8]*P*P
                + _rpc[j+2][9]*H*H          + _rpc[j+2][10]*(P*L)*H     + _rpc[j+2][11]*(L*L)*L
                + _rpc[j+2][12]*(L*P)*P     + _rpc[j+2][13]*(L*H)*H     + _rpc[j+2][14]*(L*L)*P
                + _rpc[j+2][15]*(P*P)*P     + _rpc[j+2][16]*(P*H)*H     + _rpc[j+2][17]*(L*L)*H
                + _rpc[j+2][18]*(P*P)*H     + _rpc[j+2][19]*(H*H)*H;
        }

        P       = Coeff[0]/Coeff[1];
        L       = Coeff[2]/Coeff[3];
        
        Object.m_Y = P*_rpc[1][3] + _rpc[0][3]; //Y
        Object.m_X = L*_rpc[1][2] + _rpc[0][2]; //X
        Object.m_Z = GP_H;//Z
        
        if(L < -10.0 || L > 10.0)
        {
            if(Object.m_X > 0)
                Object.m_X = Object.m_X - 360;
            else
                Object.m_X = Object.m_X + 360;
        }
        
        if(P < -10.0 || P > 10.0)
        {
            if(Object.m_Y > 0)
                Object.m_Y = Object.m_Y - 360;
            else
                Object.m_Y = Object.m_Y + 360;
            
        }
    }

    return Object;
}

void GetRayVectorFromIRPC(const double * const *IRPCs, TransParam param, const uint8 numofparam, double *imageparam, CSize imagesize, D3DPOINT &ray_vector)
{
    double minH = (double) (-1.0 * IRPCs[1][4] + IRPCs[0][4]);
    double maxH = (double) (1.0 * IRPCs[1][4] + IRPCs[0][4]);
    double midH = minH + (maxH - minH)/2.0;
    
    D2DPOINT imc(imagesize.width/2.0, imagesize.height/2.0);
    D3DPOINT s_object(0,0,midH);
    D3DPOINT e_object(0,0,-1.0 * IRPCs[1][4] + IRPCs[0][4]);
    D3DPOINT *senodes = NULL;
    
    vector<D2DPOINT> im_center;
    vector<D3DPOINT> nodes;
    
    im_center.push_back(imc);
    im_center.push_back(imc);
    nodes.push_back(s_object);
    nodes.push_back(e_object);
    
    vector<D3DPOINT> nodes_ps;
    
    senodes = GetImageHToObjectIRPC(IRPCs,2,imageparam,nodes,im_center);
    
    for(int i=0;i<nodes.size();i++)
    {
        D2DPOINT wgs_XY(senodes[i].m_X,senodes[i].m_Y);
        D2DPOINT ps_XY = wgs2ps_single(param,wgs_XY);
        
        D3DPOINT temp_ps(ps_XY.m_X,ps_XY.m_Y,senodes[i].m_Z);
        nodes_ps.push_back(temp_ps);
        
        //printf("xyz %f\t%f\t%f\n",senodes[i].m_X,senodes[i].m_Y,senodes[i].m_Z);
        //printf("xyz %f\t%f\t%f\n",temp_ps.m_X,temp_ps.m_Y,temp_ps.m_Z);
    }
    
    free(senodes);

    ray_vector.m_X = nodes_ps[1].m_X - nodes_ps[0].m_X;
    ray_vector.m_Y = nodes_ps[1].m_Y - nodes_ps[0].m_Y;
    ray_vector.m_Z = nodes_ps[1].m_Z - nodes_ps[0].m_Z;
    
    double mag = SQRT(ray_vector);
    ray_vector.m_X /= mag;
    ray_vector.m_Y /= mag;
    ray_vector.m_Z /= mag;
    //printf("ray vector IRPC %f\t%f\t%f\n",ray_vector.m_X,ray_vector.m_Y,ray_vector.m_Z);
}

void GetRayVectorFromEOBRcenter(const EO eo, CAMERA_INFO camera, CSize imagesize, double *boundary, double *minmaxH, D3DPOINT &ray_vector)
{
    double midH = minmaxH[0] + (minmaxH[1] - minmaxH[0])/2.0;
    D3DPOINT center(boundary[0] + (boundary[2] - boundary[0])/2.0, boundary[1] + (boundary[3] - boundary[1])/2.0, midH);
    
    RM M = MakeRotationMatrix(eo.m_Wl,eo.m_Pl,eo.m_Kl);
    D2DPOINT photo_coord = GetPhotoCoordinate_single(center, eo, camera, M);
    D2DPOINT image_coord = PhotoToImage_single(photo_coord, camera.m_CCDSize, imagesize);
    
    //printf("image coord %f\t%f\n",image_coord.m_X,image_coord.m_Y);
    
    ray_vector.m_X = center.m_X - eo.m_Xl;
    ray_vector.m_Y = center.m_Y - eo.m_Yl;
    ray_vector.m_Z = center.m_Z - eo.m_Zl;
    
    
    double mag = SQRT(ray_vector);
    ray_vector.m_X /= mag;
    ray_vector.m_Y /= mag;
    ray_vector.m_Z /= mag;
    
    //printf("ray vector EO center %f\t%f\t%f\n",ray_vector.m_X,ray_vector.m_Y,ray_vector.m_Z);
}

void GetAZELFromRay(const D3DPOINT ray_vector, float &AZ, float &EL)
{
    AZ = atan2(-ray_vector.m_X,-ray_vector.m_Y)*RadToDeg;
    EL = atan(fabs(ray_vector.m_Z)/sqrt(ray_vector.m_X*ray_vector.m_X + ray_vector.m_Y*ray_vector.m_Y))*RadToDeg;
    
    //printf("AZ EL %f\t%f\n",AZ, EL);
}

void GetBaseRayFromIRPC(const EO image1, const EO image2, const double * const *IRPCs1, const double * const *IRPCs2, TransParam param, const uint8 numofparam, double *imageparam, CSize imagesize1, CSize imagesize2, D3DPOINT &BR)
{
    //image1
    double maxH = image1.m_Zl;// 475000;//(double) (1.0 * IRPCs1[1][4] + IRPCs1[0][4]);
    
    D2DPOINT imc(imagesize1.width/2.0, imagesize1.height/2.0);
    D3DPOINT s_object(0,0,maxH);
    D3DPOINT *senodes = NULL;
    
    vector<D2DPOINT> im_center;
    vector<D3DPOINT> nodes;
    
    im_center.push_back(imc);
    nodes.push_back(s_object);
    
    vector<D3DPOINT> nodes_ps;
    
    senodes = GetImageHToObjectIRPC(IRPCs1,numofparam,imageparam,nodes,im_center);
    
    D2DPOINT wgs_XY(senodes[0].m_X,senodes[0].m_Y);
    D2DPOINT ps_XY = wgs2ps_single(param,wgs_XY);
    
    D3DPOINT temp_ps(ps_XY.m_X,ps_XY.m_Y,senodes[0].m_Z);
    nodes_ps.push_back(temp_ps);
    
    //printf("xyz %f\t%f\t%f\n",senodes[0].m_X,senodes[0].m_Y,senodes[0].m_Z);
    //printf("xyz %f\t%f\t%f\t%f\t%f\t%f\n",temp_ps.m_X,temp_ps.m_Y,temp_ps.m_Z,image1.m_Xl,image1.m_Yl,image1.m_Zl);
    
    free(senodes);
    
    im_center.clear();
    nodes.clear();
    
    imc.m_X = imagesize2.width/2.0;
    imc.m_Y = imagesize2.height/2.0;
    im_center.push_back(imc);
    
    maxH = image2.m_Zl;
    s_object.m_X = 0.0;
    s_object.m_Y = 0.0;
    s_object.m_Z = maxH;
    nodes.push_back(s_object);
    
    senodes = GetImageHToObjectIRPC(IRPCs2,numofparam,imageparam,nodes,im_center);
    
    wgs_XY.m_X = senodes[0].m_X;
    wgs_XY.m_Y = senodes[0].m_Y;
    ps_XY = wgs2ps_single(param,wgs_XY);
    
    temp_ps.m_X = ps_XY.m_X;
    temp_ps.m_Y = ps_XY.m_Y;
    temp_ps.m_Z = senodes[0].m_Z;
    nodes_ps.push_back(temp_ps);
    
    //printf("xyz %f\t%f\t%f\n",senodes[0].m_X,senodes[0].m_Y,senodes[0].m_Z);
    //printf("xyz %f\t%f\t%f\t%f\t%f\t%f\n",temp_ps.m_X,temp_ps.m_Y,temp_ps.m_Z,image2.m_Xl,image2.m_Yl,image2.m_Zl);
    
    free(senodes);
    
    BR.m_X = nodes_ps[1].m_X - nodes_ps[0].m_X;
    BR.m_Y = nodes_ps[1].m_Y - nodes_ps[0].m_Y;
    BR.m_Z = nodes_ps[1].m_Z - nodes_ps[0].m_Z;
    
    double mag = SQRT(BR);
    BR.m_X /= mag;
    BR.m_Y /= mag;
    BR.m_Z /= mag;
    //printf("base vector %f\t%f\t%f\n",BR.m_X,BR.m_Y,BR.m_Z);
    
}

void GetBaseRayFromIRPC(const double * const *IRPCs1, const double * const *IRPCs2, TransParam param, const uint8 numofparam, double *imageparam, CSize imagesize1, CSize imagesize2, D3DPOINT &BR)
{
    //image1
    double maxH = 475000;//(double) (1.0 * IRPCs1[1][4] + IRPCs1[0][4]);
    
    D2DPOINT imc(imagesize1.width/2.0, imagesize1.height/2.0);
    D3DPOINT s_object(0,0,maxH);
    D3DPOINT *senodes = NULL;
    
    vector<D2DPOINT> im_center;
    vector<D3DPOINT> nodes;
    
    im_center.push_back(imc);
    nodes.push_back(s_object);
    
    vector<D3DPOINT> nodes_ps;
    
    senodes = GetImageHToObjectIRPC(IRPCs1,numofparam,imageparam,nodes,im_center);
    
    D2DPOINT wgs_XY(senodes[0].m_X,senodes[0].m_Y);
    D2DPOINT ps_XY = wgs2ps_single(param,wgs_XY);
    
    D3DPOINT temp_ps(ps_XY.m_X,ps_XY.m_Y,senodes[0].m_Z);
    nodes_ps.push_back(temp_ps);
    
    //printf("xyz %f\t%f\t%f\n",senodes[0].m_X,senodes[0].m_Y,senodes[0].m_Z);
    //printf("xyz %f\t%f\t%f\t%f\t%f\t%f\n",temp_ps.m_X,temp_ps.m_Y,temp_ps.m_Z,image1.m_Xl,image1.m_Yl,image1.m_Zl);
    
    free(senodes);
    
    im_center.clear();
    nodes.clear();
    
    imc.m_X = imagesize2.width/2.0;
    imc.m_Y = imagesize2.height/2.0;
    im_center.push_back(imc);
    
    //maxH = image2.m_Zl;
    s_object.m_X = 0.0;
    s_object.m_Y = 0.0;
    s_object.m_Z = maxH;
    nodes.push_back(s_object);
    
    senodes = GetImageHToObjectIRPC(IRPCs2,numofparam,imageparam,nodes,im_center);
    
    wgs_XY.m_X = senodes[0].m_X;
    wgs_XY.m_Y = senodes[0].m_Y;
    ps_XY = wgs2ps_single(param,wgs_XY);
    
    temp_ps.m_X = ps_XY.m_X;
    temp_ps.m_Y = ps_XY.m_Y;
    temp_ps.m_Z = senodes[0].m_Z;
    nodes_ps.push_back(temp_ps);
    
    //printf("xyz %f\t%f\t%f\n",senodes[0].m_X,senodes[0].m_Y,senodes[0].m_Z);
    //printf("xyz %f\t%f\t%f\t%f\t%f\t%f\n",temp_ps.m_X,temp_ps.m_Y,temp_ps.m_Z,image2.m_Xl,image2.m_Yl,image2.m_Zl);
    
    free(senodes);
    
    BR.m_X = nodes_ps[1].m_X - nodes_ps[0].m_X;
    BR.m_Y = nodes_ps[1].m_Y - nodes_ps[0].m_Y;
    BR.m_Z = nodes_ps[1].m_Z - nodes_ps[0].m_Z;
    
    double mag = SQRT(BR);
    BR.m_X /= mag;
    BR.m_Y /= mag;
    BR.m_Z /= mag;
    //printf("base vector %f\t%f\t%f\n",BR.m_X,BR.m_Y,BR.m_Z);
    
}

void GetBaseRayFromEO(const EO image1, const EO image2, D3DPOINT &BR)
{
    BR.m_X = image2.m_Xl - image1.m_Xl;
    BR.m_Y = image2.m_Yl - image1.m_Yl;
    BR.m_Z = image2.m_Zl - image1.m_Zl;
    
    double mag = SQRT(BR);
    BR.m_X /= mag;
    BR.m_Y /= mag;
    BR.m_Z /= mag;
    //printf("base vector %f\t%f\t%f\t%f\n",BR.m_X,BR.m_Y,BR.m_Z, mag);
    
    //printf("vector 1 %f\t%f\t%f\n",ray_vector1.m_X,ray_vector1.m_Y,ray_vector1.m_Z);
    //printf("vector 2 %f\t%f\t%f\n",ray_vector2.m_X,ray_vector2.m_Y,ray_vector2.m_Z);
}

void GetPairAnglesFromRays(const double A1, const double A2, const double E1, const double E2, const D3DPOINT ray_vector1, const D3DPOINT ray_vector2, const D3DPOINT base_vector, double &CA, double &AE, double &BIE)
{
    double mag = SQRT(base_vector);
    double mag1 = SQRT(ray_vector1);
    double mag2 = SQRT(ray_vector2);
    
    double a1 = acos((ray_vector1.m_X*base_vector.m_X + ray_vector1.m_Y*base_vector.m_Y + ray_vector1.m_Z*base_vector.m_Z)/(mag*mag1))*RadToDeg;
    double a2 = 180 - acos((ray_vector2.m_X*base_vector.m_X + ray_vector2.m_Y*base_vector.m_Y + ray_vector2.m_Z*base_vector.m_Z)/(mag*mag2))*RadToDeg;
    //printf("a1 a2 CA %f\t%f\t%f\n",a1,a2,180-a1-a2);
    
    CA = fabs(180-a1-a2);
    AE = fabs((a1 - a2)/2.0);
    double BIE_num = sin(E1*DegToRad) + sin(E2*DegToRad);
    double BIE_den = sqrt(2.0)*sqrt( 1.0 + cos((A1-A2)*DegToRad) * cos(E1*DegToRad)*cos(E2*DegToRad) + sin(E1*DegToRad)*sin(E2*DegToRad) );
    BIE     = asin(BIE_num/BIE_den)*RadToDeg;
    
    //printf("vector 1 %f\t%f\t%f\n",ray_vector1.m_X,ray_vector1.m_Y,ray_vector1.m_Z);
    //printf("vector 2 %f\t%f\t%f\n",ray_vector2.m_X,ray_vector2.m_Y,ray_vector2.m_Z);
    //printf("a1 a2 CA AE BIE %f\t%f\t%f\t%f\t%f\n",a1,a2,CA,AE,BIE);
}

void GetStereoGeometryFromIRPC(const EO image1, const EO image2, const double * const *IRPCs1, const double * const *IRPCs2, const D3DPOINT ray_vector1, const D3DPOINT ray_vector2, const ImageInfo Iinfo1, const ImageInfo Iinfo2, TransParam param, const uint8 numofparam, double *imageparam, CSize imagesize1, CSize imagesize2, double &CA, double &AE, double &BIE, D3DPOINT &BR)
{
    GetBaseRayFromIRPC(image1, image2, IRPCs1, IRPCs2, param, numofparam, imageparam, imagesize1, imagesize2, BR);
    GetPairAnglesFromRays(Iinfo1.AZ_ray[0], Iinfo2.AZ_ray[0], Iinfo1.EL_ray[0], Iinfo2.EL_ray[0], ray_vector1, ray_vector2, BR, CA, AE, BIE);
}

void GetStereoGeometryFromIRPC(const double * const *IRPCs1, const double * const *IRPCs2, const D3DPOINT ray_vector1, const D3DPOINT ray_vector2, const ImageInfo Iinfo1, const ImageInfo Iinfo2, TransParam param, const uint8 numofparam, double *imageparam, CSize imagesize1, CSize imagesize2, double &CA, double &AE, double &BIE, D3DPOINT &BR)
{
    GetBaseRayFromIRPC(IRPCs1, IRPCs2, param, numofparam, imageparam, imagesize1, imagesize2, BR);
    GetPairAnglesFromRays(Iinfo1.AZ_ray[0], Iinfo2.AZ_ray[0], Iinfo1.EL_ray[0], Iinfo2.EL_ray[0], ray_vector1, ray_vector2, BR, CA, AE, BIE);
}

void GetStereoGeometryFromEO(const EO image1, const EO image2, const D3DPOINT ray_vector1, const D3DPOINT ray_vector2, const ImageInfo Iinfo1, const ImageInfo Iinfo2, double &CA, double &AE, double &BIE, D3DPOINT &BR)
{
    GetBaseRayFromEO(image1, image2, BR);
    GetPairAnglesFromRays(Iinfo1.AZ_ray[1], Iinfo2.AZ_ray[1], Iinfo1.EL_ray[1], Iinfo2.EL_ray[1], ray_vector1, ray_vector2, BR, CA, AE, BIE);
}

D2DPOINT GetObjectToImageRPC_single_mpp(const double * const *_rpc, const uint8 _numofparam, const double *_imageparam, D3DPOINT _GP)
{
    double L, P, H, Line, Samp, deltaP, deltaR;
    double Coeff[4];

    L       = (_GP.m_X - _rpc[0][2])/_rpc[1][2];
    P       = (_GP.m_Y - _rpc[0][3])/_rpc[1][3];
    H       = (_GP.m_Z - _rpc[0][4])/_rpc[1][4];

    if(L < -10.0 || L > 10.0)
    {
        if(_GP.m_X > 0)
            _GP.m_X = _GP.m_X - 360;
        else
            _GP.m_X = _GP.m_X + 360;

        L       = (_GP.m_X - _rpc[0][2])/_rpc[1][2];
    }

    if(P < -10.0 || P > 10.0)
    {
        if(_GP.m_Y > 0)
            _GP.m_Y = _GP.m_Y - 360;
        else
            _GP.m_Y = _GP.m_Y + 360;

        P       = (_GP.m_Y - _rpc[0][3])/_rpc[1][3];
    }

    for(int j=0;j<4;j++)
    {
        Coeff[j]    = _rpc[j+2][0]*1.0          + _rpc[j+2][1]*L            + _rpc[j+2][2]*P
            + _rpc[j+2][3]*H            + _rpc[j+2][4]*L*P          + _rpc[j+2][5]*L*H
            + _rpc[j+2][6]*P*H          + _rpc[j+2][7]*L*L          + _rpc[j+2][8]*P*P
            + _rpc[j+2][9]*H*H          + _rpc[j+2][10]*(P*L)*H     + _rpc[j+2][11]*(L*L)*L
            + _rpc[j+2][12]*(L*P)*P     + _rpc[j+2][13]*(L*H)*H     + _rpc[j+2][14]*(L*L)*P
            + _rpc[j+2][15]*(P*P)*P     + _rpc[j+2][16]*(P*H)*H     + _rpc[j+2][17]*(L*L)*H
            + _rpc[j+2][18]*(P*P)*H     + _rpc[j+2][19]*(H*H)*H;
    }

    Line     = ((Coeff[0]/Coeff[1])*_rpc[1][0] + _rpc[0][0]); //Line
    Samp     = ((Coeff[2]/Coeff[3])*_rpc[1][1] + _rpc[0][1]); //Sample

    deltaP      = _imageparam[0];
    deltaR      = _imageparam[1];

    //printf("deltaP R %f\t%f\n",deltaP,deltaR);
    D2DPOINT IP;
    IP.m_Y      = deltaP + Line;
    IP.m_X      = deltaR + Samp;

    return IP;
}

D2DPOINT* GetObjectToImage(uint16 _numofpts, D2DPOINT *_GP, double *boundary, double imageres)
{
    D2DPOINT *IP;
    
    IP        = (D2DPOINT*)malloc(sizeof(D2DPOINT)*_numofpts);
    
#pragma omp parallel for schedule(guided)
    for(int i=0;i<_numofpts;i++)
    {
        IP[i].m_X = ( _GP[i].m_X - boundary[0])/imageres;
        IP[i].m_Y = (-_GP[i].m_Y + boundary[3])/imageres;
    }
    
    return IP;
}

D2DPOINT GetObjectToImage_single(uint16 _numofpts, D2DPOINT _GP, double *boundary, double imageres)
{
    D2DPOINT IP;
    
    IP.m_X = ( _GP.m_X - boundary[0])/imageres;
    IP.m_Y = (-_GP.m_Y + boundary[3])/imageres;
    
    return IP;
}



D2DPOINT* OriginalToPyramid(uint16 numofpts,D2DPOINT* InCoord, D2DPOINT Startpos, uint8 Pyramid_step)
{
    int i;
    D2DPOINT* out;

    out = (D2DPOINT*)malloc(sizeof(D2DPOINT)*numofpts);
    for(i=0;i<numofpts;i++)
    {
        out[i].m_X      = (InCoord[i].m_X/pwrtwo(Pyramid_step)) - Startpos.m_X;
        out[i].m_Y      = (InCoord[i].m_Y/pwrtwo(Pyramid_step)) - Startpos.m_Y;
    }

    return out;
    
}

D2DPOINT* PyramidToOriginal(uint16 numofpts,D2DPOINT* InCoord, D2DPOINT Startpos, uint8 Pyramid_step)
{
    int i;
    D2DPOINT* out;
    out = (D2DPOINT*)malloc(sizeof(D2DPOINT)*numofpts);

    if(Pyramid_step > 0)
    {
        for(i=0;i<numofpts;i++)
        {
            out[i].m_X      = (InCoord[i].m_X +  Startpos.m_X)*pwrtwo(Pyramid_step) + pwrtwo(Pyramid_step)/2.0;
            out[i].m_Y      = (InCoord[i].m_Y +  Startpos.m_Y)*pwrtwo(Pyramid_step) + pwrtwo(Pyramid_step)/2.0;
        }
    }
    else
    {
        for(i=0;i<numofpts;i++)
        {
            out[i].m_X      = InCoord[i].m_X +  Startpos.m_X;
            out[i].m_Y      = InCoord[i].m_Y +  Startpos.m_Y;
        }
    }

    return out;
}


//Aerial Frame Camera
RM MakeRotationMatrix(double o, double p, double k)
{
    double M[9] = {0.0};
    RM m_rm;
    o *= DegToRad;
    p *= DegToRad;
    k *= DegToRad;
    
    m_rm.m11 = cos(p)*cos(k);
    m_rm.m12 = sin(o)*sin(p)*cos(k)+cos(o)*sin(k);
    m_rm.m13 = -cos(o)*sin(p)*cos(k)+sin(o)*sin(k);
    
    m_rm.m21 = -cos(p)*sin(k);
    m_rm.m22 = -sin(o)*sin(p)*sin(k)+cos(o)*cos(k);
    m_rm.m23 = cos(o)*sin(p)*sin(k)+sin(o)*cos(k);
    
    m_rm.m31 = sin(p);
    m_rm.m32 = -sin(o)*cos(p);
    m_rm.m33 = cos(o)*cos(p);
    
    return m_rm;
}

D2DPOINT *GetPhotoCoordinate(D3DPOINT *A, EO Photo, int _numofpts, CAMERA_INFO Camera, RM M)
{
    D2DPOINT *IP;
    
    IP      = (D2DPOINT*)malloc(sizeof(D2DPOINT)*_numofpts);
    
#pragma omp parallel for schedule(guided)
    for(int i=0;i<_numofpts;i++)
    {
        double q = M.m31*(A[i].m_X - Photo.m_Xl) + M.m32*(A[i].m_Y - Photo.m_Yl) + M.m33*(A[i].m_Z - Photo.m_Zl);
        double r = M.m11*(A[i].m_X - Photo.m_Xl) + M.m12*(A[i].m_Y - Photo.m_Yl) + M.m13*(A[i].m_Z - Photo.m_Zl);
        double s = M.m21*(A[i].m_X - Photo.m_Xl) + M.m22*(A[i].m_Y - Photo.m_Yl) + M.m23*(A[i].m_Z - Photo.m_Zl);
        
        IP[i].m_X = Camera.m_ppx - Camera.m_focalLength*(r/q);
        IP[i].m_Y = Camera.m_ppy - Camera.m_focalLength*(s/q);
        
        double RR;
        double x_, y_;
        x_ = IP[i].m_X - Camera.m_ppx;
        y_ = IP[i].m_Y - Camera.m_ppy;
        RR = (x_*x_)+(y_*y_);
        
        double dxR, dxT;
        double dyR, dyT;
        double dR_R;
        dR_R = Camera.k1*RR + Camera.k2*RR*RR + Camera.k3*RR*RR*RR;
        dxR = x_*dR_R;
        dyR = y_*dR_R;
        dxT = Camera.p1*(RR + 2*x_*x_) + 2*Camera.p2*x_*y_;
        dyT = Camera.p2*(RR + 2*y_*y_) + 2*Camera.p1*x_*y_;
        
        IP[i].m_X = IP[i].m_X - dxR - dxT;
        IP[i].m_Y = IP[i].m_Y - dyR - dyT - Camera.a1*x_ - Camera.a2*y_;
        
        /*
        double q = M.m31*(A[i].m_X - Photo.m_Xl) + M.m32*(A[i].m_Y - Photo.m_Yl) + M.m33*(A[i].m_Z - Photo.m_Zl);
        double r = M.m11*(A[i].m_X - Photo.m_Xl) + M.m12*(A[i].m_Y - Photo.m_Yl) + M.m13*(A[i].m_Z - Photo.m_Zl);
        double s = M.m21*(A[i].m_X - Photo.m_Xl) + M.m22*(A[i].m_Y - Photo.m_Yl) + M.m23*(A[i].m_Z - Photo.m_Zl);
        
        IP[i].m_X = Camera.m_ppx - Camera.m_focalLength*(r/q);
        IP[i].m_Y = Camera.m_ppy - Camera.m_focalLength*(s/q);
         */
    }
    
    return IP;
}

D3DPOINT *GetObjectCoordinate(D2DPOINT *a, double z,EO Photo, int _numofpts, CAMERA_INFO Camera, RM M)
{
    D3DPOINT *A;
    A      = (D3DPOINT*)malloc(sizeof(D3DPOINT)*_numofpts);
    
#pragma omp parallel for schedule(guided)
    for(int i=0;i<_numofpts;i++)
    {
        D2DPOINT corr_pt;
        
        double RR;
        double x_, y_;
        x_ = a[i].m_X - Camera.m_ppx;
        y_ = a[i].m_Y - Camera.m_ppy;
        RR = (x_*x_)+(y_*y_);
        
        double dxR, dxT;
        double dyR, dyT;
        double dR_R;
        dR_R = Camera.k1*RR + Camera.k2*RR*RR + Camera.k3*RR*RR*RR;
        dxR = x_*dR_R;
        dyR = y_*dR_R;
        dxT = Camera.p1*(RR + 2*x_*x_) + 2*Camera.p2*x_*y_;
        dyT = Camera.p2*(RR + 2*y_*y_) + 2*Camera.p1*x_*y_;
        
        corr_pt.m_X = x_ + dxR + dxT;
        corr_pt.m_Y = y_ + dyR + dyT + Camera.a1*x_ + Camera.a2*y_;
        
        A[i].m_Z = z;
        A[i].m_X =  (A[i].m_Z - Photo.m_Zl)*(M.m11*(corr_pt.m_X) + M.m21*(corr_pt.m_Y) + M.m31*(-Camera.m_focalLength))
        /(M.m13*(corr_pt.m_X) + M.m23*(corr_pt.m_Y) + M.m33*(-Camera.m_focalLength)) + Photo.m_Xl;
        A[i].m_Y =  (A[i].m_Z - Photo.m_Zl)*(M.m12*(corr_pt.m_X) + M.m22*(corr_pt.m_Y) + M.m32*(-Camera.m_focalLength))
        /(M.m13*(corr_pt.m_X) + M.m23*(corr_pt.m_Y) + M.m33*(-Camera.m_focalLength)) + Photo.m_Yl;
        
        
        /*
        A[i].m_Z = z;
        A[i].m_X =  (A[i].m_Z - Photo.m_Zl)*(M.m11*(a[i].m_X - Camera.m_ppx) + M.m21*(a[i].m_Y - Camera.m_ppy) + M.m31*(-Camera.m_focalLength))
        /(M.m13*(a[i].m_X - Camera.m_ppx) + M.m23*(a[i].m_Y - Camera.m_ppy) + M.m33*(-Camera.m_focalLength)) + Photo.m_Xl;
        A[i].m_Y =  (A[i].m_Z - Photo.m_Zl)*(M.m12*(a[i].m_X - Camera.m_ppx) + M.m22*(a[i].m_Y - Camera.m_ppy) + M.m32*(-Camera.m_focalLength))
        /(M.m13*(a[i].m_X - Camera.m_ppx) + M.m23*(a[i].m_Y - Camera.m_ppy) + M.m33*(-Camera.m_focalLength)) + Photo.m_Yl;
         */
    }
    
    return A;
}

D2DPOINT *PhotoToImage(D2DPOINT *_photo, int _numofpts, float _CCDSize, CSize _imgsize)
{
    D2DPOINT *m_ImageCoord;
    m_ImageCoord      = (D2DPOINT*)malloc(sizeof(D2DPOINT)*_numofpts);
    
#pragma omp parallel for schedule(guided)
    for(int i=0;i<_numofpts;i++)
    {
        m_ImageCoord[i].m_X = ( _photo[i].m_X + _imgsize.width*(_CCDSize*UMToMM)/2.)/(_CCDSize*UMToMM);
        m_ImageCoord[i].m_Y = (-_photo[i].m_Y + _imgsize.height*(_CCDSize*UMToMM)/2.)/(_CCDSize*UMToMM);
    }
    
    return m_ImageCoord;
}

D2DPOINT *ImageToPhoto(D2DPOINT *_image, int _numofpts, float _CCDSize, CSize _imgsize)
{
    D2DPOINT *m_PhotoCoord;
    m_PhotoCoord      = (D2DPOINT*)malloc(sizeof(D2DPOINT)*_numofpts);
    
#pragma omp parallel for schedule(guided)
    for(int i=0;i<_numofpts;i++)
    {
        m_PhotoCoord[i].m_X =   _image[i].m_X*(_CCDSize*UMToMM) - _imgsize.width*(_CCDSize*UMToMM)/2.;
        m_PhotoCoord[i].m_Y = -(_image[i].m_Y*(_CCDSize*UMToMM) - _imgsize.height*(_CCDSize*UMToMM)/2.);
    }
    
    return m_PhotoCoord;
}

D2DPOINT GetPhotoCoordinate_single(const D3DPOINT A, const EO Photo, const CAMERA_INFO Camera, const RM M)
{
    D2DPOINT IP;
    
    double q = M.m31*(A.m_X - Photo.m_Xl) + M.m32*(A.m_Y - Photo.m_Yl) + M.m33*(A.m_Z - Photo.m_Zl);
    double r = M.m11*(A.m_X - Photo.m_Xl) + M.m12*(A.m_Y - Photo.m_Yl) + M.m13*(A.m_Z - Photo.m_Zl);
    double s = M.m21*(A.m_X - Photo.m_Xl) + M.m22*(A.m_Y - Photo.m_Yl) + M.m23*(A.m_Z - Photo.m_Zl);
    
    IP.m_X = Camera.m_ppx - Camera.m_focalLength*(r/q);
    IP.m_Y = Camera.m_ppy - Camera.m_focalLength*(s/q);
    
    double RR;
    double x_, y_;
    x_ = IP.m_X - Camera.m_ppx;
    y_ = IP.m_Y - Camera.m_ppy;
    RR = (x_*x_)+(y_*y_);
    
    double dxR, dxT;
    double dyR, dyT;
    double dR_R;
    dR_R = Camera.k1*RR + Camera.k2*RR*RR + Camera.k3*RR*RR*RR;
    dxR = x_*dR_R;
    dyR = y_*dR_R;
    dxT = Camera.p1*(RR + 2*x_*x_) + 2*Camera.p2*x_*y_;
    dyT = Camera.p2*(RR + 2*y_*y_) + 2*Camera.p1*x_*y_;
    
    IP.m_X = IP.m_X - dxR - dxT;
    IP.m_Y = IP.m_Y - dyR - dyT - Camera.a1*x_ - Camera.a2*y_;
    
    return IP;
}

D3DPOINT GetObjectCoordinate_single(D2DPOINT a, double z,EO Photo, CAMERA_INFO Camera, RM M)
{
    D3DPOINT A;
    D2DPOINT corr_pt;
    
    double RR;
    double x_, y_;
    x_ = a.m_X - Camera.m_ppx;
    y_ = a.m_Y - Camera.m_ppy;
    RR = (x_*x_)+(y_*y_);
    
    double dxR, dxT;
    double dyR, dyT;
    double dR_R;
    dR_R = Camera.k1*RR + Camera.k2*RR*RR + Camera.k3*RR*RR*RR;
    dxR = x_*dR_R;
    dyR = y_*dR_R;
    dxT = Camera.p1*(RR + 2*x_*x_) + 2*Camera.p2*x_*y_;
    dyT = Camera.p2*(RR + 2*y_*y_) + 2*Camera.p1*x_*y_;
    
    corr_pt.m_X = x_ + dxR + dxT;
    corr_pt.m_Y = y_ + dyR + dyT + Camera.a1*x_ + Camera.a2*y_;
    
    A.m_Z = z;
    A.m_X =  (A.m_Z - Photo.m_Zl)*(M.m11*(corr_pt.m_X) + M.m21*(corr_pt.m_Y) + M.m31*(-Camera.m_focalLength))
    /(M.m13*(corr_pt.m_X) + M.m23*(corr_pt.m_Y) + M.m33*(-Camera.m_focalLength)) + Photo.m_Xl;
    A.m_Y =  (A.m_Z - Photo.m_Zl)*(M.m12*(corr_pt.m_X) + M.m22*(corr_pt.m_Y) + M.m32*(-Camera.m_focalLength))
    /(M.m13*(corr_pt.m_X) + M.m23*(corr_pt.m_Y) + M.m33*(-Camera.m_focalLength)) + Photo.m_Yl;
    
    /*
    A.m_X =  (A.m_Z - Photo.m_Zl)*(M.m11*(a.m_X - Camera.m_ppx) + M.m21*(a.m_Y - Camera.m_ppy) + M.m31*(-Camera.m_focalLength))
    /(M.m13*(a.m_X - Camera.m_ppx) + M.m23*(a.m_Y - Camera.m_ppy) + M.m33*(-Camera.m_focalLength)) + Photo.m_Xl;
    A.m_Y =  (A.m_Z - Photo.m_Zl)*(M.m12*(a.m_X - Camera.m_ppx) + M.m22*(a.m_Y - Camera.m_ppy) + M.m31*(-Camera.m_focalLength))
    /(M.m13*(a.m_X - Camera.m_ppx) + M.m23*(a.m_Y - Camera.m_ppy) + M.m33*(-Camera.m_focalLength)) + Photo.m_Yl;
     */
    return A;
}

D2DPOINT PhotoToImage_single(const D2DPOINT _photo, const double _CCDSize, const CSize _imgsize)
{
    D2DPOINT m_ImageCoord;
    
    m_ImageCoord.m_X = ( _photo.m_X + _imgsize.width*(_CCDSize*UMToMM)/2.)/(_CCDSize*UMToMM);
    m_ImageCoord.m_Y = (-_photo.m_Y + _imgsize.height*(_CCDSize*UMToMM)/2.)/(_CCDSize*UMToMM);
    
    return m_ImageCoord;
}

D2DPOINT ImageToPhoto_single(D2DPOINT _image, float _CCDSize, CSize _imgsize)
{
    D2DPOINT m_PhotoCoord;
    
    m_PhotoCoord.m_X =   _image.m_X*(_CCDSize*UMToMM) - _imgsize.width*(_CCDSize*UMToMM)/2.;
    m_PhotoCoord.m_Y = -(_image.m_Y*(_CCDSize*UMToMM) - _imgsize.height*(_CCDSize*UMToMM)/2.);
    
    return m_PhotoCoord;
}
