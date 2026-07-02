void COls_calcDoc::OnFOUR_ADD() 
{
//フーリエ係数の集計
	// TODO: この位置にコマンド ハンドラ用のコードを追加してください
	FCOMP cell;
	FCOMP cell2;
	FCOMP cellx1;
	FCOMP cellx2;

	double *FFT1_r;
	double *FFT1_i;
	double *pDFT1_r;
	double *pDFT1_i;
	double *pDFT2_r;
	double *pDFT2_i;
	double *Ang;
	double *Radius;
	double *Rp_SG;//Savitzky Golay 2次多項式 平滑化微分の計算結果が入る
	int SG;//平滑化微分点数=2xSG+1
	SG=5;
	double Rp_border;//ﾋﾟｰｸ検知のしきい値
	Rp_border=(double)0.2;//ここでは生data MAX値の1/5に設定

	double *Rp_frq;//検知されたﾋﾟｰｸの周波数が入る
	double *Rp_spec;//検知されたﾋﾟｰｸのｽﾍﾟｸﾄﾙが入る
	int Rp_num;//ﾋﾟｰｸの個数
	Rp_num=0;
	(double *)FFT1_r=new double[fft_size*2];
	(double *)FFT1_i=new double[fft_size*2];
	(double *)pDFT1_r=new double[fft_size*fft_size];
	(double *)pDFT1_i=new double[fft_size*fft_size];
	(double *)pDFT2_r=new double[fft_size*fft_size];
	(double *)pDFT2_i=new double[fft_size*fft_size];
	(double *)Ang=new double[181];
	(double *)Radius=new double[fft_size/2+1];
	(double *)Rp_SG=new double[fft_size/2];
	(double *)Rp_frq=new double[fft_size/4];
	(double *)Rp_spec=new double[fft_size/4];
	int i,j,k;
	for(i=0;i<181;i++){
		Ang[i]=0;
	}
	for(i=0;i<fft_size/2+1;i++){
		Radius[i]=0;
	}
//----------FFT演算部--------------------------

	int ret;//bit反転の結果が入る変数
	int m;//0 or 1計算結果を代入する変数の区別
	int n;//0 or 1計算される変数の区別
	int l;//FFT stageの変数
	int bat;//ﾊﾞﾀﾌﾗｲ演算の単位
	int sec;//ﾊﾞﾀﾌﾗｲ演算の回数
	for(i=0;i<fft_size;i++){//行
		for(k=0;k<fft_size/2;k++){
//			FFT1_r[0*fft_size+k]=(float)*(LSzData+fft_size*i+k)+(float)*(LSzData+fft_size*i+k+fft_size/2);
//			FFT1_r[0*fft_size+k]=(double)*(LSzDATA+PixelsX*(i+offy)+k+offx)+(double)*(LSzDATA+PixelsX*(i+offy)+(k+offx)+fft_size/2);
			FFT1_r[0 * fft_size + k] = (double)*(LSzDATA + fft_size*i + k) + (double)*(LSzDATA + fft_size*i + k + fft_size / 2);
			FFT1_i[0*fft_size+k]=0;
//			cellx2=(FCOMP)polar((float)1,(float)(-3.1415926535*2/fft_size*k))*((float)*(LSzData+fft_size*i+k)-(float)*(LSzData+fft_size*i+k+fft_size/2));
//			cellx2=(FCOMP)polar((double)1,(double)(-3.1415926535*2/fft_size*k))*((double)*(LSzDATA+PixelsX*(i+offy)+k+offx)-(double)*(LSzDATA+PixelsX*(i+offy)+(k+offx)+fft_size/2));
			cellx2 = (FCOMP)polar((double)1, (double)(-3.1415926535 * 2 / fft_size*k))*((double)*(LSzDATA + fft_size*i + k ) - (double)*(LSzDATA + fft_size*i + k  + fft_size / 2));
			FFT1_r[0*fft_size+k+fft_size/2]=cellx2.real();
			FFT1_i[0*fft_size+k+fft_size/2]=cellx2.imag();
		}
		for(l=1;l<fft_stage;l++){
			m=(int)l%2;
			n=(int)((m+1)%2);
			bat=(int)(fft_size/pow((double)2,(double)l)/2);
			sec=(int)pow((double)2,(double)l);
			for(j=0;j<sec;j++){
				for(k=0;k<bat;k++){
					cell=FCOMP((double)FFT1_r[n*fft_size+j*bat*2+k],(double)FFT1_i[n*fft_size+j*bat*2+k]);
					cell2=FCOMP((double)FFT1_r[n*fft_size+j*bat*2+bat+k],(double)FFT1_i[n*fft_size+j*bat*2+bat+k]);
					cellx1=cell+cell2;
					FFT1_r[m*fft_size+j*bat*2+k]=cellx1.real();
					FFT1_i[m*fft_size+j*bat*2+k]=cellx1.imag();
					cellx2=(FCOMP)polar((double)1,(double)(-3.1415926535/bat*k))*((cell)-(cell2));
					FFT1_r[m*fft_size+j*bat*2+bat+k]=cellx2.real();
					FFT1_i[m*fft_size+j*bat*2+bat+k]=cellx2.imag();
				}
			}
		}
		for(j=0;j<fft_size;j++){//列
			ret=bit_return(j,fft_stage);
			pDFT1_r[i*fft_size+ret]=FFT1_r[m*fft_size+j];
			pDFT1_i[i*fft_size+ret]=FFT1_i[m*fft_size+j];
		}
	}//ここまで終了98.5.21

	for(i=0;i<fft_size;i++){//列
		for(k=0;k<fft_size/2;k++){
			cell=FCOMP((double)pDFT1_r[k*fft_size+i],(double)pDFT1_i[k*fft_size+i]);
			cell2=FCOMP((double)pDFT1_r[(k+fft_size/2)*fft_size+i],(double)pDFT1_i[(k+fft_size/2)*fft_size+i]);
			cellx1=(cell)+(cell2);
			FFT1_r[0*fft_size+k]=cellx1.real();
			FFT1_i[0*fft_size+k]=cellx1.imag();
			cellx2=(FCOMP)polar((double)1,(double)(-3.1415926535*2/fft_size*k))*((cell)-(cell2));
			FFT1_r[0*fft_size+k+fft_size/2]=cellx2.real();
			FFT1_i[0*fft_size+k+fft_size/2]=cellx2.imag();
		}
		for(l=1;l<fft_stage;l++){
			m=(int)l%2;
			n=(int)((m+1)%2);
			bat=(int)(fft_size/pow((double)2,(double)l)/2);
			sec=(int)pow((double)2,(double)l);
			for(j=0;j<sec;j++){
				for(k=0;k<bat;k++){
					cell=FCOMP((double)FFT1_r[n*fft_size+j*bat*2+k],(double)FFT1_i[n*fft_size+j*bat*2+k]);
					cell2=FCOMP((double)FFT1_r[n*fft_size+j*bat*2+bat+k],(double)FFT1_i[n*fft_size+j*bat*2+bat+k]);
					cellx1=(cell)+(cell2);
					FFT1_r[m*fft_size+j*bat*2+k]=cellx1.real();
					FFT1_i[m*fft_size+j*bat*2+k]=cellx1.imag();
					cellx2=(FCOMP)polar((double)1,(double)(-3.1415926535/bat*k))*((cell)-(cell2));
					FFT1_r[m*fft_size+j*bat*2+bat+k]=cellx2.real();
					FFT1_i[m*fft_size+j*bat*2+bat+k]=cellx2.imag();
				}
			}
		}
		for(j=0;j<fft_size;j++){//列
			ret=bit_return(j,fft_stage);
			pDFT2_r[ret*fft_size+i]=FFT1_r[m*fft_size+j];
			pDFT2_i[ret*fft_size+i]=FFT1_i[m*fft_size+j];
		}
	}
//--------ここからsylk作成-----------------
	CString slk;
	char slk2[200000];
	char str2[41];
	char str3[41];
	char str4[41];
	strcpy(slk2,"");
	int Rows;
	if(fft_size/2<180){
		Rows=180;
	}else{
		Rows=fft_size/2;
	}
	Rows+=120;
	_itoa(Rows,str2,10);
	
	slk="ID;PWXL\nB;Y";
	slk+=str2;
	slk+=";X10\nC;Y1;X1;K\"";
	slk+=FILE_NAME;
	slk+="\"\nC;Y2;X1;K\"";
	slk+=Com1;
	slk+="\"\n";
	slk+="C;Y3;X1;K\"";
	slk+=Com2;
	slk+="\"\n";
	slk+="C;Y4;X1;K\"";
	slk+=s_date;
	slk+="\"\n";
	slk+="C;Y5;X1;K\"";
	slk+=s_time;
	slk+="\"\n";

	slk+="C;Y7;X1;K\"XYcalib\"\n";
	slk+="C;Y8;X1;K\"Zcoeff\"\n";
	_gcvt((double)XYcalib,8,str2);
	slk+="C;Y7;X2;K";
	slk+=str2;
	slk+='\n';
	_gcvt((double)Zcoeff,8,str2);
	slk+="C;Y8;X2;K";
	slk+=str2;
	slk+='\n';

	slk+="C;Y9;X1;K\"XSize\"\n";
	slk+="C;Y10;X1;K\"YSize\"\n";
	_gcvt((double)dSize.cx,8,str2);
	slk+="C;Y9;X2;K";
	slk+=str2;
	slk+='\n';
	_gcvt((double)dSize.cy,8,str2);
	slk+="C;Y10;X2;K";
	slk+=str2;
	slk+='\n';

	slk+="C;Y11;X1;K\"XYunit\"\n";
	slk+="C;Y12;X1;K\"Zunit\"\n";
	slk+="C;Y11;X2;K\"";
	slk+=XYunit;
	slk+="\"\n";
	slk+="C;Y12;X2;K\"";
	slk+=Zunit;
	slk+="\"\n";

	slk+="C;Y14;X1;K\"MR1\"\n";
	slk+="C;Y15;X1;K\"MR2\"\n";
	_gcvt((double)O_dat.MR1,8,str2);
	slk+="C;Y14;X2;K";
	slk+=str2;
	slk+='\n';
	_gcvt((double)O_dat.MR2,8,str2);
	slk+="C;Y15;X2;K";
	slk+=str2;
	slk+='\n';
	slk+="C;Y16;X1;K\"Rpk\"\n";
	slk+="C;Y17;X1;K\"Rk\"\n";
	slk+="C;Y18;X1;K\"Rvk\"\n";
	_gcvt((double)O_dat.Rpk,8,str2);
	slk+="C;Y16;X2;K";
	slk+=str2;
	slk+='\n';
	_gcvt((double)O_dat.Rk,8,str2);
	slk+="C;Y17;X2;K";
	slk+=str2;
	slk+='\n';
	_gcvt((double)O_dat.Rvk,8,str2);
	slk+="C;Y18;X2;K";
	slk+=str2;
	slk+='\n';

	slk += "C;Y20;X1;K\"DFT size\"\n";
	slk += "C;Y20;X2;K\"";
	_itoa(fft_size, str2, 10);
	slk += str2;
	slk += " x ";
	slk += str2;
	slk += "\"\n";


//-----ここから集計ｽﾀｰﾄ--------------
	double spec;
	int x,x2;
	int y,y2;
	int r;
	int ang;
	for(i=0;i<fft_size-1;i++){//i:X
		x=(fft_size/2+1+i)%fft_size;
		x2=(fft_size/2+1+i)-fft_size;
		for(j=0;j<fft_size/2-1;j++){//j:Y
			y=(fft_size/2+1+j)%fft_size;
			y2=(fft_size/2+1+j)-fft_size;
			r=(int)(pow((pow((double)x2,(double)2)+pow((double)y2,(double)2)),0.5)+0.5);
			if(r<fft_size/2 && r>0){
				ang=(int)((double)(acos(x2/pow((x2*x2+y2*y2),0.5))/3.141592*180)+0.5);
//0702ここ殺した				spec=(double)(Zcoeff*pow(pow((pDFT2_r[y*fft_size+x]),2)+pow((pDFT2_i[y*fft_size+x]),2),0.5)/fft_size/fft_size);
				spec=(double)(pow(pow((pDFT2_r[y*fft_size+x]),2)+pow((pDFT2_i[y*fft_size+x]),2),0.5)/fft_size/fft_size);	//0702 Zcoeffを削除
				Ang[ang]+=spec*2;
				Radius[r]+=spec*2;
			}else if(r==0){
//0702ここ殺した				spec=(double)(Zcoeff*pow(pow((pDFT2_r[y*fft_size+x]),2)+pow((pDFT2_i[y*fft_size+x]),2),0.5)/fft_size/fft_size);
				spec=(double)(pow(pow((pDFT2_r[y*fft_size+x]),2)+pow((pDFT2_i[y*fft_size+x]),2),0.5)/fft_size/fft_size);	//0702 Zcoeffを削除
//				Radius[r]=spec;
				Radius[r]=0;
			}
		}
//-------------y=0の処理---------------
		j=fft_size/2-1;
		y=(fft_size/2+1+j)%fft_size;
		y2=(fft_size/2+1+j)-fft_size;
		r=(int)(pow((pow((double)x2,(double)2)+pow((double)y2,(double)2)),0.5)+0.5);
		if(r<fft_size/2 && r>0){
			ang=(int)((double)(acos(x2/pow((x2*x2+y2*y2),0.5))/3.141592*180)+0.5);
//0702殺し			spec=(double)(Zcoeff*pow(pow((pDFT2_r[y*fft_size+x]),2)+pow((pDFT2_i[y*fft_size+x]),2),0.5)/fft_size/fft_size);
			spec=(double)(pow(pow((pDFT2_r[y*fft_size+x]),2)+pow((pDFT2_i[y*fft_size+x]),2),0.5)/fft_size/fft_size);//0702 Zcoeff削除
			Ang[ang]+=spec;
			Radius[r]+=spec;
		}else if(r==0){
//0702殺し			spec=(double)(Zcoeff*pow(pow((pDFT2_r[y*fft_size+x]),2)+pow((pDFT2_i[y*fft_size+x]),2),0.5)/fft_size/fft_size);
			spec=(double)(pow(pow((pDFT2_r[y*fft_size+x]),2)+pow((pDFT2_i[y*fft_size+x]),2),0.5)/fft_size/fft_size);//0702 Zcoeff排除
//			Radius[r]=spec;
			Radius[r]=0;
		}
//-------------y=0の処理終了-----------
		
	}
	Ang[0]+=Ang[180];
	Ang[180]=Ang[0];
//-------集計終了----------------------------
//-------ﾋﾟｰｸ算出ｽﾀｰﾄ------------------------
//Savitzky Golayの2次多項式 平滑化微分
	double Rpmax;
	double Rp3;//ﾋﾟｰｸ高さ(3点平均)
	Rpmax=0;
	for(i=1;i<fft_size/2;i++){
		if(Radius[i]>Rpmax){Rpmax=Radius[i];}
	}
	Rpmax=Rpmax*Rp_border;//しきい値の設定
	double bf;//前点の微分係数
	bf=-1;//最初の点の評価のため負にｾｯﾄ
	double haba;//ﾋﾟｰｸ点幅
	double frq;
	haba=0;
	//-------平滑化できない低周波の波形でのﾋﾟｰｸの算出---------
	//生波形で山となっており、かつしきい値より大きいものはすべてﾋﾟｰｸとして数える
	for(i=1;i<SG+1;i++){
		if((Radius[i-1]<Radius[i])&&(Radius[i+1]<Radius[i])){
			if(Radius[i]>Rpmax){
				Rp_frq[Rp_num]=(double)i;
				Rp_spec[Rp_num]=Radius[i];
				Rp_num+=1;
			}
		}
	}
	//--------------------------------------------------------
	for(i=SG+1;i<fft_size/2-SG;i++){
			Rp_SG[i]=0;
		for(j=1;j<=SG;j++){
			Rp_SG[i]+=j*(Radius[i+j]-Radius[i-j]);
		}
		if(Rp_SG[i]==0){
			if(bf<0){
				haba=0;
				continue;
			}else{
				haba+=1;
				continue;
			}
		}
		if(Rp_SG[i]>0){haba=0;}
		if(Rp_SG[i]<0){
			if(bf==0){j=i-1-(int)(haba/2);frq=i-1-(int)(haba/2);}
			if(bf>0){j=i-(int)(haba/2+0.5);frq=i-haba-1+((haba)/2)+(bf/(bf-Rp_SG[i]));}
			if(bf<0){
				haba=0;
				bf=Rp_SG[i];
				continue;
			}
			//-----ﾋﾟｰｸ高さ計算 しきい値以上ならﾋﾟｰｸとして検出する------
			Rp3=(Radius[j-1]+Radius[j]+Radius[j+1])/3;
			if(Rp3>Rpmax){
				Rp_frq[Rp_num]=frq;
				Rp_spec[Rp_num]=Rp3;
				Rp_num+=1;
			}
		}
		bf=Rp_SG[i];
	}

//-------ﾋﾟｰｸ算出終了------------------------
	for(i=0;i<fft_size/2;i++){
		_itoa(1,str3,10);
		_itoa(i,str4,10);
		_itoa(i+22,str2,10);
		strcat(slk2,"C;Y");
		strcat(slk2,str2);
		strcat(slk2,";X");
		strcat(slk2,str3);
		strcat(slk2,";K");
		strcat(slk2,str4);
		strcat(slk2,"\n");
		_itoa(2,str3,10);
		_gcvt(Radius[i],8,str4);
		strcat(slk2,"C;Y");
		strcat(slk2,str2);
		strcat(slk2,";X");
		strcat(slk2,str3);
		strcat(slk2,";K");
		strcat(slk2,str4);
		strcat(slk2,"\n");
	}
	for(i=0;i<=180;i++){
		_itoa(6,str3,10);
		_itoa(i,str4,10);
		_itoa(i+22,str2,10);
		strcat(slk2,"C;Y");
		strcat(slk2,str2);
		strcat(slk2,";X");
		strcat(slk2,str3);
		strcat(slk2,";K");
		strcat(slk2,str4);
		strcat(slk2,"\n");
		_itoa(7,str3,10);
		_gcvt(Ang[i],8,str4);
		strcat(slk2,"C;Y");
		strcat(slk2,str2);
		strcat(slk2,";X");
		strcat(slk2,str3);
		strcat(slk2,";K");
		strcat(slk2,str4);
		strcat(slk2,"\n");
	}
//	int le=strlen(slk2);
	slk+=slk2;
//--------ﾋﾟｰｸ値のSYLK化-------------------------------
		_itoa(fft_size/2+24,str2,10);
		strcpy(slk2,"");
		strcat(slk2,"C;Y");
		strcat(slk2,str2);
		strcat(slk2,";X1;K\"周波数ﾋﾟｰｸ解析\"\n");
		_itoa(fft_size/2+25,str2,10);
		strcat(slk2,"C;Y");
		strcat(slk2,str2);
		strcat(slk2,";X1;K\"演算方法\"\n");
		strcat(slk2,"C;Y");
		strcat(slk2,str2);
		strcat(slk2,";X2;K\"平滑化微分\"\n");
		_itoa(2*SG+1,str3,10);
		strcat(slk2,"C;Y");
		strcat(slk2,str2);
		strcat(slk2,";X3;K\"");
		strcat(slk2,str3);
		strcat(slk2,"点\"\n");
		_itoa(fft_size/2+26,str2,10);
		strcat(slk2,"C;Y");
		strcat(slk2,str2);
		strcat(slk2,";X1;K\"ﾋﾟｰｸ高さ\"\n");
		strcat(slk2,"C;Y");
		strcat(slk2,str2);
		strcat(slk2,";X2;K\"単純平均\"\n");
		_itoa(3,str3,10);
		strcat(slk2,"C;Y");
		strcat(slk2,str2);
		strcat(slk2,";X3;K\"");
		strcat(slk2,str3);
		strcat(slk2,"点\"\n");
		_itoa(fft_size/2+27,str2,10);
		strcat(slk2,"C;Y");
		strcat(slk2,str2);
		strcat(slk2,";X1;KNo.\n");
		strcat(slk2,"C;Y");
		strcat(slk2,str2);
		strcat(slk2,";X2;K\"次数\"\n");
		strcat(slk2,"C;Y");
		strcat(slk2,str2);
		strcat(slk2,";X3;K\"波長\"\n");
		strcat(slk2,"C;Y");
		strcat(slk2,str2);
		strcat(slk2,";X4;K\"高さ\"\n");
		strcat(slk2,"C;Y");
		strcat(slk2,str2);
		strcat(slk2,";X5;K\"平滑化\"\n");
		for(i=0;i<Rp_num;i++){
			_itoa(fft_size/2+27+i+1,str2,10);
			_itoa(i+1,str3,10);
			strcat(slk2,"C;Y");
			strcat(slk2,str2);
			strcat(slk2,";X1;K");
			strcat(slk2,str3);
			strcat(slk2,"\n");
//			_itoa(Rp_frq[i],str3,10);
			_gcvt(Rp_frq[i],8,str3);
			strcat(slk2,"C;Y");
			strcat(slk2,str2);
			strcat(slk2,";X2;K");
			strcat(slk2,str3);
			strcat(slk2,"\n");
			_gcvt((XYcalib*fft_size/Rp_frq[i]),8,str3);
			strcat(slk2,"C;Y");
			strcat(slk2,str2);
			strcat(slk2,";X3;K");
			strcat(slk2,str3);
			strcat(slk2,"\n");
			_gcvt(Rp_spec[i],8,str3);
			strcat(slk2,"C;Y");
			strcat(slk2,str2);
			strcat(slk2,";X4;K");
			strcat(slk2,str3);
			strcat(slk2,"\n");
			if(Rp_frq[i]<=SG+1){
				strcat(slk2,"C;Y");
				strcat(slk2,str2);
				strcat(slk2,";X5;K\"無し\"\n");
			}
		}
		delete[] FFT1_r;
		delete[] FFT1_i;
		delete[] pDFT1_r;
		delete[] pDFT1_i;
		delete[] pDFT2_r;
		delete[] pDFT2_i;
		delete[] Ang;
		delete[] Radius;

		delete[] Rp_SG;
		delete[] Rp_frq;
		delete[] Rp_spec;

//--------ﾋﾟｰｸ値のSYLK化終了---------------------------
	slk+=slk2;
	slk+="E\n";
		int len=strlen(slk.GetBuffer(10));
		int len2;
		len2=len+1;
		HGLOBAL hGmem;
			void *lpGmem;
			UINT flg;
			flg=(GMEM_MOVEABLE|GMEM_ZEROINIT);
			hGmem=GlobalAlloc(flg,DWORD(len2));
			if(hGmem==NULL){TRACE("error=%d\n",GetLastError());}
			lpGmem=GlobalLock(hGmem);
			memcpy(lpGmem,slk.GetBuffer(10),len);
			GlobalUnlock(hGmem);
		if (OpenClipboard(NULL))
		{
			BeginWaitCursor();
			EmptyClipboard();
			SetClipboardData (CF_SYLK, hGmem );
			CloseClipboard();
			EndWaitCursor();
		}

//	GlobalFree(hGmem);

}
