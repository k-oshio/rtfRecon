//
// Filt project (image filter)
//


#import <RecKit/RecKit.h>
#import "RecImagePW.h"
#import <RecKit/timer_macros.h>

#import <NumKit/NumKit.h>

void	test0();	// logistic
void	test1();	// Pinwheel (vs DLR)
void	test2();	// noise estimation from 256 ave (NCI/se-epi/070919/raw01/IMG_base)
void	test3();	// gaussian vs. lorentzian (ratio)
void	test4();	// DLR analysis (difference)

void		pinwheelFilter(RecImage *coef, float sd, float flr, int mode);
// 1: multiplicative, 2: additive
void		filtChar1WithCoef(RecImage *hist, RecImage *coefi, RecImage *coefo, int dim, float gain);
void		filtChar2WithCoef(RecImage *hist, RecImage *coefi, RecImage *coefo, int dim, float gain);

void		pinwheelCoefHist(RecImage *hist, RecImage *coef, float gain, float max, BOOL lg);
RecImage	*magVar(RecImage *coef);	// measure of noise correlation

int
main()
{
    @autoreleasepool {
	//	test0();
	//	test1();
		test2();
	//	test3();
	//	test4();
    }
	return 0;
}

void
test0()
{
	/* COVID-19
report# china	global	global delta
	58	81116	
	57	81116	
	56	81077	86434	13874
	55	81048	72469	10955
	54	81021	61518	9751
	53	80991	51767	7488
	52	80981	44067	6703
	51	80955	37371	4596
	50	80924	32778	4105
	49	80904	28673	3948
	48	80859	24727	3610
	47	80813	21110	3633
	46	80711	17481	2727
	45	80565	14768	2089
	44	80422	12668	2013
	43	80304	10566	1792
	42	80174	8774	1600
	41	79968	7169	1160
	40	79394	6009	1318
	39	78961	4691	1027
	38	78630	3664	746
	37	78191	2918	459
	36	77780	2459	390
	35	77262	2069	300
	34	77042	1769	367
	33	76392	1402	202
	32	75569	1200	127
	31	74675	1073	149
	30	74280	924		120
	29	72528	804		10
	28	70635	794		111
	27	51174	683		157
	26	50054	526		21
	25	48548	505		58
	24	46550	447		6
	23	44730	441		46
	22	42708	395		76
	21	40235	319		12
	20	37251	307		19
	19	34598	288		18
	18	31211	270		54
	17	28060	216		25
	16	24363	191		32
	15	20471	159		6
	14	17238	153		7
	13	14471	146		14
	12	11821	132		26
	*/
	/* COVID-19
report# italy	delta	usa		delta
	72	105702	4053	163199	22559
	71	101739	4050	140640	17987
	70	97689	5217	122653	19332
	69	92472	5974	103321	18093
	68	86498	5959	85228	16894
	67	80539	6153	68334	4764
	66	74386	5210	63570	11656
	65	69176	5249	51914	9750
	64	63927	4789	42164	10591
	63	59138	5560	31573	16354
	62	53578	6557	15219	0
	61	47021	5986	15219	4777
	60	41035	5322	10442	3335
	59	35713	4207	7078	3551
	58	31506	3526	3536	1822
	57	27980	3233	3503	1825
	*/
	/*
report# korea	delta
	72	9887	101	
	71	9786	125	
	70	9661	78
	69	9583	105
	68	9478	146
	67	9332	91
	66	9241	104
	65	9137	100
	64	9037	76
	63	8961	64
	62	8897	98
	61	8799	147
	60	8652	239
	59	8413	58
	58	8355	35
	57	8320	84
	56	8236	74
	55	8162	76
	54	8086	107
	53	7979	110
	52	7869	114
	51	7755	242
	50	7513	131
	49	7382	248
	48	7134	367
	47	6767	483
	46	6284	518
	45	5766	438
	44	5328	516
	43	4812	600
	42	4212	476
	41	3736	586
	40	3150	813
	39	2337	571
	38	1766	505
	37	1261	284
	36	977		214
	35	763		161
	34	602		256
	33	346		142
	32	204		13
	31	104		53
	30	51		20
	*/
	/* Influenza -2019
	2020	959.4400	23
	2019	959.4400	57
	2018	838.0034	55
	2017	735.5818	39
	2016	637.3564	40
	2015	499.3689	39
	2014	240.9158	34
	2013	131.4437	37
	2012	142.5100	43
	2011	104.3246	32
	2010	141.2875	40
	2009	100.085		37
	2008	100.0416	18
	2007	 94.2439	33
	2006	 81.1675	32
	2005	 65.2820	50
	
	*/
	/* Influenza 2019-2020 (x10000)
	13		92
	12		71
	11		75
	10		73
	9		81
	8		100
	7		105
	6		88
	5		75
	4		72
	3		69
	2		72
	1		75
	*/


	int		i;
	float	x, y1, y2, y3, y4;
	float	dy1, dy2, dy3, dy4;
	float	mx1 = 81000;	// china
	float	ctr1 = 19;	// 2/7/2020
	float	r1 = 4.5;
	float	mx2 = 80000;	// korea	
	float	ctr2 = 51;	// 3/11/2020
	float	r2 = 4.5;
	float	mx3 = 125000;	// italy
	float	ctr3 = 64;	// 3/24/2020
	float	r3 = 4.9;
	float	mx4 = 500000;	// USA
	float	ctr4 = 75;	// 4/5/2020
	float	r4 = 4.5;
	
	for (i = 10; i < 80; i++) {
		x = i;
		y1 = 1.0 / (1 + exp(- (x - ctr1) / r1));
		y2 = 1.0 / (1 + exp(- (x - ctr2) / r2));
		y3 = 1.0 / (1 + exp(- (x - ctr3) / r3));
		y4 = 1.0 / (1 + exp(- (x - ctr4) / r4));
		dy1 = y1 * (1.0 - y1);
		dy2 = y2 * (1.0 - y2);
		dy3 = y3 * (1.0 - y3);
		dy4 = y4 * (1.0 - y4);
	//	printf("%5.2f %5.2f %5.2f\n", x, mx1 * y1, mx2 * y2);
	//	printf("%5.2f %5.2f %5.2f %5.2f %5.2f\n", x, mx1 * y1, mx2 * y2, mx1 / r1 * dy1, mx2 / r2 * dy2);
	//	printf("%5.2f %5.2f %5.2f %5.2f\n", x, mx1 * y1, mx2 * y2, mx2 / r2 * dy2);
		printf("%5.2f %5.2f %5.2f %5.2f %5.2f\n", x, mx3 * y3, mx3 / r3 * dy3, mx4 * y4, mx4 / r4 * dy4);
	}
}

void
test1()
{
	RecImage	*img = [RecImage imageWithKOImage:@"CANON_DLR_Sub1/img.img"];
	RecImage	*slc1, *slc2;
	RecImage	*coef1, *coef2, *dif;
	RecImage	*filt, *fSlc;
	RecImage	*hist, *hSlc;
	RecImage	*correl;
	RecLoop		*zLp, *xLp;
	int			i, n = 5;
	int			order = 6;
	int			dim = 1024;
	RecImage	*kl, *kern;
	float		sd = 0.002, mx;

	system("rm img_*");

// pinwheel kernel
	kl = Rec_pinwheel_param(order);
	kern = Rec_pinwheel_kernel(img, kl);	// freq domain
	zLp = [RecLoop loopWithDataLength:n];
	xLp = [RecLoop loopWithDataLength:dim];

// DL filter
	slc1 = [img sliceAtIndex:0];	// 1 AVE
	coef1 = Rec_pinwheel_decomp(slc1, kern);
	mx = [coef1 maxVal];

	fSlc = [slc1 copy];
	hSlc = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:[coef1 zLoop], xLp, nil];

	pinwheelCoefHist(hSlc, coef1, 10.0, mx, YES);
	[hSlc saveAsKOImage:@"img_hist0"];

	hist = [hSlc copy];
	[hist addLoop:zLp];
	filt = [fSlc copy];
	[filt addLoop:zLp];
	correl = [filt copy];
	for (i = 0; i < n; i++) {
		slc2 = [img sliceAtIndex:i + 1];
		coef2 = Rec_pinwheel_decomp(slc2, kern);
	//	filtChar1WithCoef(hSlc, coef1, coef2, dim, 10.0);
		filtChar2WithCoef(hSlc, coef1, coef2, dim, 10.0);	// ???? ####
		[hist copySlice:hSlc atIndex:i];
		dif = [coef2 copy];
		[dif subImage:coef1];
	//	dif = [dif sdForLoop:[dif zLoop]];
		dif = magVar(dif);

		[correl copySlice:dif atIndex:i];
	}
	[hist saveAsKOImage:@"img_hist1"];
	[img saveAsKOImage:@"img_filt1"];
	[correl saveAsKOImage:@"img_corr1"];	// not done yet

// pinwheel filter
	for (i = 0; i < n; i++) {
		coef2 = [coef1 copy];
		pinwheelFilter(coef2, sd * (i + 1), 0.9, 2);
	//	pinwheelFilter(coef2, sd * (i + 1), 1.0, 3);	// 1: hard threshold, 2: gaussian, 3:lorentzian

		dif = [coef2 copy];
		[dif subImage:coef1];	// coeff set, take var from here
	//	dif = [dif sdForLoop:[dif zLoop]];
		dif = magVar(dif);
		[correl copySlice:dif atIndex:i];

	//	filtChar1WithCoef(hSlc, coef1, coef2, dim, 10.0);
		filtChar2WithCoef(hSlc, coef1, coef2, dim, 10.0);
		[hist copySlice:hSlc atIndex:i];
		
		fSlc = Rec_pinwheel_synth(coef2, kern);
		[filt copySlice:fSlc atIndex:i];
	}
	[hist saveAsKOImage:@"img_hist2"];
	[filt saveAsKOImage:@"img_filt2"];
	[correl saveAsKOImage:@"img_corr2"];	// not done yet
}

// ###
// *1. select reference (~ 100 avg)
// 2. calc MSE w.r.t. ref
void
test2()
{
	RecImage	*img = [RecImage imageWithKOImage:@"SE_EPI/img_base"];
	RecImage	*mask;
	RecImage	*img_p;
	RecImage	*avg, *slc, *sum, *tmp;
	RecImage	*ref, *dif, *filt, *correl;
//	NSString	*path;
	RecLoop		*zLp, *xLp;
	int			i, len, n = 100;
	float		rms, rms0;
	RecImage	*kl, *kern;
	RecImage	*coef1, *coef2;
	int			order = 6;
	float		sd = 0.04, sd0 = sd/8.0;
	float		gain = 0.6;
	int			flt = 2;		// 0:pass through, 1:hard thres, 2:gaussian, 3:lorentzian
	int			mode = 1;

// pinwheel filter
	kl = Rec_pinwheel_param(order);
	kern = Rec_pinwheel_kernel(img, kl);	// freq domain
	Rec_pinwheel_dump_kernel(order);

// acquired image
	switch (mode) {
	case 0:
		len = [img zDim];
		sum = [img sliceAtIndex:0];
		avg = [img copy];
		[avg crop:[avg zLoop] to:n];
		for (i = 1; i < len; i++) {
			[sum addImage:[img sliceAtIndex:i]];
			slc = [sum copy];
			slc = [slc multByConst:1.0/(i + 1)];
			[avg copySlice:slc atIndex:i];
		}
		[avg saveAsKOImage:@"img_ave"];

		dif = [avg copy];
		ref = [avg sliceAtIndex:n - 1];
		[dif subImage:ref];
		[dif saveAsKOImage:@"img_dif"];
		rms0 = [[dif sliceAtIndex:0] rmsVal];
		for (i = 0; i < n; i++) {
			tmp = [dif sliceAtIndex:i];
			printf("%f %f\n", sqrt((float)i), [tmp rmsVal]/rms0);
		}
		break;

// generated image
	case 1:
		len = [img zDim];
		ref = [img avgForLoop:[img zLoop]];
		[ref saveAsKOImage:@"img_ref"];
		
		img_p = [RecImage imageWithImage:img];
		[img_p copyImage:ref];
		[img_p addGWN:sd relative:YES];
	//	[img_p magnitude];
		[img_p saveAsKOImage:@"img_phantom"];
		img = [img_p copy];
		avg = [img copy];
		[avg crop:[avg zLoop] to:n];
		sum = [img sliceAtIndex:0];

		for (i = 1; i < n; i++) {
			[sum addImage:[img sliceAtIndex:i]];
			slc = [sum copy];
			slc = [slc multByConst:1.0/(i + 1)];
			[avg copySlice:slc atIndex:i];
		}
//	[avg magnitude];
		[avg saveAsKOImage:@"img_avg"];
		slc = [avg sliceAtIndex:n-1];
		coef1 = Rec_pinwheel_decomp(slc, kern);
		slc = [avg sliceAtIndex:0];
		coef2 = Rec_pinwheel_decomp(slc, kern);
		[coef2 subImage:coef1];
		[coef2 saveAsKOImage:@"img_noise"];

		xLp = [RecLoop loopWithDataLength:512];
		slc = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:[coef1 zLoop], xLp, nil];
		pinwheelCoefHist(slc, coef2, 1.0, [slc maxVal], NO);
		[slc saveAsKOImage:@"img_nz_hist"];
		


		dif = [avg copy];
		[dif subImage:ref];
		[dif saveAsKOImage:@"img_dif"];
		rms0 = [[dif sliceAtIndex:0] rmsVal];
		for (i = 0; i < n; i++) {
			tmp = [dif sliceAtIndex:i];
			printf("%f %f\n", sqrt((float)i), [tmp rmsVal]/rms0);
		}
		break;
	}
	mask = [ref copy];
	[mask thresAt:0.05];
	[mask saveAsKOImage:@"img_mask"];
//exit(0);

	slc = [dif sliceAtIndex:0];
	[slc multByImage:mask];
	rms0 = [slc rmsVal];	// rms val of added noise

	slc = [img sliceAtIndex:1];
//	[slc magnitude];	// phantom has neg noise
	coef1 = Rec_pinwheel_decomp(slc, kern);

	zLp = [RecLoop loopWithDataLength:10];
	filt = [RecImage imageOfType:RECIMAGE_REAL withLoops:zLp, [coef1 yLoop], [coef1 xLoop], nil];
	dif = [filt copy];
	correl = [filt copy];
	for (i = 0; i < 10; i++) {
		coef2 = [coef1 copy];
		sd = sd0 * (i + 1);
		pinwheelFilter(coef2, sd, gain, flt);	// 1: hard threshold, 2: gaussian, 3:lorentzian
		slc = Rec_pinwheel_synth(coef2, kern);
		[filt copySlice:slc atIndex:i];
		[slc subImage:ref];
		[dif copySlice:slc atIndex:i];
		[slc multByImage:mask];
		rms = [slc rmsVal] / rms0;
	//	printf("%d %f %f", i, sd, rms);
		printf("%d %f", i, rms);

	// noise correlation measure
		[coef2 subImage:coef1];
	[correl copySlice:coef2 atIndex:i];
		rms = [coef2 maxVal];
		[coef2 multByConst:1.0/rms];
		rms = [coef2 rmsVal];
		printf(" %f\n", rms);
	}
	[filt saveAsKOImage:@"img_filt"];
	[dif saveAsKOImage:@"img_mse"];
	[correl saveAsKOImage:@"img_corr"];
}

void
pinwheelFilter(RecImage *coef, float sd, float mix, int mode)
{
	int			i, j;
	int			imgSize, len;
	float		*pp, val, mx, fr;
	float		fw;			// freq / magnitude weight
	RecImage	*gfactor;	// coil gain
	float		x, y, r, *gp;
	int			xDim, yDim;

	gfactor = [RecImage imageOfType:RECIMAGE_REAL withLoops:[coef yLoop], [coef xLoop], nil];
	xDim = [gfactor xDim];
	yDim = [gfactor yDim];
	gp = [gfactor data];
	for (i = 0; i < yDim; i++) {
		y = 2.0 * ((float)i - yDim/2) / yDim;
		for (j = 0; j < xDim; j++) {
			x = 2.0 * ((float)j - xDim/2) / xDim;
			r = x*x + y*y;
			gp[i * xDim + j] = 4.0 * exp(- r / 1.0);
		}
	}

// PW filter
	pp = [coef real];	// real only
	imgSize = [coef xDim] * [coef yDim];
	len = [coef zDim];
	mx = [coef maxVal];
	for (i = 0; i < imgSize; i++) {
		for (j = 0; j < len; j++) {
			val = pp[j * imgSize + i];
			fr = fabs(val/mx);
			switch (mode) {
			case 0 :	// pass-through
				break;
			case 1 :	// hard threshold
				if (fr < sd) {
					val *= 1.0 - mix;
				}
				break;
			case 2 :	// gaussian + floor
				fw = exp(-fr * fr / (2 * sd * sd));
				fw *= mix;
				fw = 1.0 - fw;
				val *= fw;
				break;
			case 3 :	// wiener ### not done yet
				fw = sd*sd / (fr*fr + sd*sd);
				fw *= mix;
				fw = 1.0 - fw;
				val *= fw;
				break;
			}
			pp[j * imgSize + i] = val;
		}
	}
}

// filter characteristics... av of coef multiplier
void
filtChar1WithCoef(RecImage *hist, RecImage *coefi, RecImage *coefo, int dim, float gain)
{
	RecImage	*cnt;
	RecLoop		*xLp;
	int			i, ixi, ixh, k, n;
	int			zDim;
	int			x, y, sign;
	float		mx;
	float		*pi, *po, *pc, *pp; // *qq;
	float		vo, vi;
	
	xLp = [RecLoop loopWithDataLength:dim];

	[hist clear];
	cnt = [hist copy];
	zDim = [coefi zDim];
	pp = [hist real];
//	qq = [hist imag];
	pc = [cnt data];
	mx = [coefi maxVal];

// 2D filter param diagram
	pi = [coefi data];
	po = [coefo data];
	n = [coefi xDim] * [coefi yDim];

// multiplier histo (re:avg, im:var)
	for (k = 1; k < zDim; k++) {	// skip DC
		y = k;
		for (i = 0; i < n; i++) {
			ixi = k * n + i;	// image index
		//	val = po[ix];
			x = dim/2 + 0.5 * gain * pi[ixi] * dim / mx;
			ixh = y * dim + x;	// hist index
			if (x >= 0 && x < dim) {
				vi = pi[ixi];
				vo = po[ixi];
				if (vi > 0) {
					sign = 1;
				} else {
					sign = -1;
				}
				pp[ixh] += fabs(vo);
			//	qq[ixh] += (vo - vi) * sign; // * (vi - vo);
				pc[ixh] += fabs(vi);
			}
		}
	}
	for (i = 0; i < dim * zDim; i++) {
		if (pc[i] > 0) {
			pp[i] /= pc[i];
		//	qq[i] /= pc[i];
		//} else {
		//	pp[i] = 1.0;
		}
	}
}

// filter characteristics 2... av / var of coef difference
void
filtChar2WithCoef(RecImage *hist, RecImage *coefi, RecImage *coefo, int dim, float gain)
{
	RecImage	*cnt;
	RecLoop		*xLp;
	int			i, ixi, ixh, k, n;
	int			zDim;
	int			x, y, sign;
	float		mx;
	float		*pi, *po, *pc, *pp, *qq;
	float		vo, vi;

	xLp = [RecLoop loopWithDataLength:dim];

	[hist clear];
	cnt = [hist copy];
	zDim = [coefi zDim];
	pp = [hist real];
	qq = [hist imag];
	pc = [cnt data];
	mx = [coefi maxVal];

// 2D filter param diagram
	pi = [coefi data];
	po = [coefo data];
	n = [coefi xDim] * [coefi yDim];

// difference histo == mean
	for (k = 1; k < zDim; k++) {	// skip DC
		y = k;
		for (i = 0; i < n; i++) {
			ixi = k * n + i;	// image index
			x = dim/2 + 0.5 * gain * pi[ixi] * dim / mx;
			ixh = y * dim + x;	// hist index
			if (x >= 0 && x < dim) {
				vi = pi[ixi];
				vo = po[ixi];
				if (vi > 0) {
					sign = 1;
				} else {
					sign = -1;
				}
			//	pp[ixh] += (vo - vi) * sign;
				pp[ixh] += (fabs(vo) - fabs(vi));
				pc[ixh] += fabs(vi);
			}
		}
	}
// difference histo == var
	for (k = 1; k < zDim; k++) {	// skip DC
		y = k;
		for (i = 0; i < n; i++) {
			ixi = k * n + i;	// image index
			x = dim/2 + 0.5 * gain * pi[ixi] * dim / mx;
			ixh = y * dim + x;	// hist index
			if (x >= 0 && x < dim) {
				vi = pi[ixi];
				vo = po[ixi];
				vo = fabs(vo) - fabs(vi) - pp[ixh];
				qq[ixh] += vo * vo;	// ## not correct
			}
		}
	}
	for (i = 0; i < dim * zDim; i++) {
		if (pc[i] > 0) {
			pp[i] *= 100.0 / pc[i];
			qq[i] *= 1.0 / pc[i];
			qq[i] = sqrt(qq[i]);
		}
	}
}

void
pinwheelCoefHist(RecImage *hist, RecImage *coef, float gain, float max, BOOL lg)
{
	int			i, ix, k, n;
	int			zDim, dim;
	int			x, y;
	float		val;
	float		mx;
	float		*p, *pp;

	[hist clear];

	zDim = [coef zDim];
	dim = [hist xDim];
	pp = [hist data];
	if (max == 0) {
		mx = [coef maxVal];
	} else {
		mx = max;
	}

// 2D filter param diagram
	p = [coef data];
	n = [coef xDim] * [coef yDim];
	for (k = 1; k < zDim; k++) {	// skip DC
		y = k;
		for (i = 0; i < n; i++) {
			ix = k * n + i;
			val = p[ix];
			x = dim/2 - 0.5 * gain * p[ix] * dim / mx;
			if (x < 0) x = 0;
			if (x >= dim) x = dim-1;
			pp[y * dim + x] += 1.0;
		}
	}
	for (x = 0; x < dim; x++) {
		for (y = 1; y < zDim; y++) {
			pp[x] += pp[y * dim + x] / zDim;
		}
	}
	if (lg) {
		[hist logP1];
	}
}

RecImage *
magVar(RecImage *coef)	// measure of noise correlation
{
	RecImage	*var;
	RecImage	*tmp = [coef copy];
	int			i, j, n = [coef xDim] * [coef yDim];
	int			zDim = [coef zDim];
	float		*p, *pp, mx, mn, val;

	var = [RecImage imageOfType:RECIMAGE_REAL withLoops:[coef yLoop], [coef xLoop], nil];
	[tmp magnitude];
	
	p = [tmp data];
	pp = [var data];

	for (i = 0; i < n; i++) {
		mx = mn = p[i];
		for (j = 1; j < zDim; j++) {
			val = p[i + j * n];
			if (mx < val) {
				mx = val;
			} else
			if (val < mn) {
				mn = val;
			}
		}
		pp[i] = mx - mn;
	}
	return var;
}

void
test3()
{
	int		i, n = 64;
	float	x, y1, y2, sd = 5;

	for (i = 0; i < n; i++) {
		x = (float)i - n/2;
		y1 = exp(- x * x / (2 * sd * sd));	// gaussian
		y2 = sd * sd / (x * x + sd * sd);	// lorentzian
		printf("%f %f %f\n", x, y1, y2);
	}
}

void
test4()
{
	RecImage	*img = [RecImage imageWithKOImage:@"CANON_DLR_Sub1/img.img"];
	RecImage	*slc1, *slc2;
	RecImage	*coef1, *coef2, *dif;
	RecImage	*filt, *fSlc;
	RecImage	*hist, *hSlc;
	RecImage	*correl;
	RecLoop		*zLp, *xLp;
	int			i, n = 5;
	int			order = 6;
	int			dim = 1024;
	RecImage	*kl, *kern;
	float		mx;

	system("rm img_*");

// pinwheel kernel
	kl = Rec_pinwheel_param(order);
	kern = Rec_pinwheel_kernel(img, kl);	// freq domain
	zLp = [RecLoop loopWithDataLength:n];
	xLp = [RecLoop loopWithDataLength:dim];

// DL filter
	slc1 = [img sliceAtIndex:0];	// 1 AVE
	coef1 = Rec_pinwheel_decomp(slc1, kern);
	mx = [coef1 maxVal];

	fSlc = [slc1 copy];
	hSlc = [RecImage imageOfType:RECIMAGE_REAL withLoops:[coef1 zLoop], xLp, nil];

	pinwheelCoefHist(hSlc, coef1, 10.0, mx, YES);
	[hSlc saveAsKOImage:@"img_hist0"];

	hist = [hSlc copy];
	[hist addLoop:zLp];
	filt = [fSlc copy];
	[filt addLoop:zLp];
	correl = [filt copy];

	for (i = 0; i < n; i++) {
		slc2 = [img sliceAtIndex:i + 1];
		[slc2 subImage:slc1];
		[img copySlice:slc2 atIndex:i + 1];
		coef2 = Rec_pinwheel_decomp(slc2, kern);
	//	filtChar1WithCoef(hSlc, coef1, coef2, dim, 10.0);
		pinwheelCoefHist(hSlc, coef2, 10.0, mx, YES);
		[hist copySlice:hSlc atIndex:i];
		dif = [coef2 copy];
		[dif subImage:coef1];
	//	dif = [dif sdForLoop:[dif zLoop]];
		dif = magVar(dif);

		[correl copySlice:dif atIndex:i];
	}
	[hist saveAsKOImage:@"img_hist1"];
	[img saveAsKOImage:@"img_diff"];
	[correl saveAsKOImage:@"img_corr1"];	// not done yet

}