#include "../INCLUDE/globalMatrix.h"

#include <../INCLUDE/Mjday.h>
#include <../INCLUDE/SAT_Const.h>
#include <string.h>

Matrix Cnm(181,181);
Matrix Snm(181,181);
Matrix obs(46,4);

void GGM03S() {
    FILE *fid = fopen("../DATA/GGM03S.txt", "r");
    if (fid==NULL) {
        printf("Error opening file\n");
        exit(-1);
    }
    double aux;
    for (int i=1;i<=181;i++) {
        for (int j=1;j<=i;j++) {
            fscanf(fid,"%lf %lf %lf %lf %lf %lf", &aux, &aux,&Cnm(i,j),&Snm(i,j),&aux,&aux);
        }
    }
    fclose(fid);
}

void GEOS3() {
    FILE *fid = fopen("../DATA/GEOS3.txt", "r");
    if (fid==NULL) {
        printf("Error opening file\n");
        exit(-1);
    }
    int Y,MO,D,H,MI,S;
    double AZ,EL,DIST;
    char tline[57],y[5],mo[3],d[3],h[3],mi[3],s[6],az[10],el[9],dist[10],aux[2];
    for (int i=1; i<=46;i++) {
        fgets(tline, sizeof(tline),fid);
        strncpy(y,&tline[0],4);
        y[4]='\0';
        Y=atoi(y);

        strncpy(mo,&(tline[5]),2);
        mo[2]='\0';
        MO=atoi(mo);

        strncpy(d,&(tline[8]),2);
        d[2]='\0';
        D=atoi(d);

        strncpy(h,&(tline[12]),2);
        h[2]='\0';
        H=atoi(h);

        strncpy(mi,&(tline[15]),2);
        mi[2]='\0';
        MI=atoi(mi);

        strncpy(s,&(tline[18]),5);
        s[5]='\0';
        S=atof(s);

        strncpy(az,&(tline[25]),9);
        az[9]='\0';
        AZ=atof(az);

        strncpy(el,&(tline[35]),8);
        el[8]='\0';
        EL=atof(el);

        strncpy(dist,&(tline[44]),9);
        dist[9]='\0';
        DIST=atof(dist);

        obs(i,1) = Mjday(Y,MO,D,H,MI,S);
        obs(i,2) = consts.Rad*AZ;
        obs(i,3) = consts.Rad*EL;
        obs(i,4) = 1e3*DIST;
    }
    fclose(fid);
}
