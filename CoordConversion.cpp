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
        A[i].m_Z = z;
        A[i].m_X =  (A[i].m_Z - Photo.m_Zl)*(M.m11*(a[i].m_X - Camera.m_ppx) + M.m21*(a[i].m_Y - Camera.m_ppy) + M.m31*(-Camera.m_focalLength))
        /(M.m13*(a[i].m_X - Camera.m_ppx) + M.m23*(a[i].m_Y - Camera.m_ppy) + M.m33*(-Camera.m_focalLength)) + Photo.m_Xl;
        A[i].m_Y =  (A[i].m_Z - Photo.m_Zl)*(M.m12*(a[i].m_X - Camera.m_ppx) + M.m22*(a[i].m_Y - Camera.m_ppy) + M.m31*(-Camera.m_focalLength))
        /(M.m13*(a[i].m_X - Camera.m_ppx) + M.m23*(a[i].m_Y - Camera.m_ppy) + M.m33*(-Camera.m_focalLength)) + Photo.m_Yl;
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

    
    return IP;
}

D3DPOINT GetObjectCoordinate_single(D2DPOINT a, double z,EO Photo, CAMERA_INFO Camera, RM M)
{
    D3DPOINT A;
    
    A.m_Z = z;
    A.m_X =  (A.m_Z - Photo.m_Zl)*(M.m11*(a.m_X - Camera.m_ppx) + M.m21*(a.m_Y - Camera.m_ppy) + M.m31*(-Camera.m_focalLength))
    /(M.m13*(a.m_X - Camera.m_ppx) + M.m23*(a.m_Y - Camera.m_ppy) + M.m33*(-Camera.m_focalLength)) + Photo.m_Xl;
    A.m_Y =  (A.m_Z - Photo.m_Zl)*(M.m12*(a.m_X - Camera.m_ppx) + M.m22*(a.m_Y - Camera.m_ppy) + M.m31*(-Camera.m_focalLength))
    /(M.m13*(a.m_X - Camera.m_ppx) + M.m23*(a.m_Y - Camera.m_ppy) + M.m33*(-Camera.m_focalLength)) + Photo.m_Yl;

    
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
