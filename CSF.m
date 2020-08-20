//
// CSF project
//


#import <RecKit/RecKit.h>
#import "RecImagePW.h"
#import <RecKit/timer_macros.h>

#import <NumKit/NumKit.h>

// util
void		mgphs_to_cpx1(NSString *base, NSString *inser, NSString *ser, int nslc);	// no b0
void		mgphs_to_cpx2(NSString *base, NSString *inser, NSString *ser, int nslc);	// b0 / bxxx

// pre-processing (mag/phs -> cpx, single slice)
void	rinkan_1();
void	rinkan_2();
void	rinkan_3();
void	rinkan_4();
void	rinkan_5();
void	nasu_1();
void	read_csv();
void	reslice();
void	read_raw_aqp();
void	read_pc();

// mg PCA
void	mg_pca();
// phs PCA
void	phs_pca();
// combined PCA
void	cpx_pca();

// mag/phs freq
void	mp_freq();

// cardiac phase avg
void	card_avg();
void	card_avg2();	// rewrite ...
void	calc_div();		// divergence (compression)

// IVIM (T2) sim
void	t2sim();

// PC flow analysis
void	pc_roi();
void	pc_avg();

// IVCM (intra-voxel coherent motion)
void	ivcm();

// quantitative time-SLIP
void    ts_1();
void    t1sim();    // time-SLIP param est
void    partial();  // partial volume effect for T2 est

int
main()
{
    @autoreleasepool {
//		rinkan_1();
//		rinkan_4();
//		rinkan_5();
//		mg_pca();
//		phs_pca();
//		nasu_1();
//		mp_freq();
//		card_avg();
//		card_avg2();
//		cpx_pca();
//		read_csv();
//		reslice();
//		read_raw_aqp();
//		read_pc();
//		pc_roi();
		pc_avg();
//		t2sim();
//		calc_div();
//		ivcm();
//        ts_1();
//        t1sim();
//        partial();
    }
	return 0;
}

float
phs_scale(RecImage *phs)
{
	float	mx, mn;
	mx = [phs maxVal];
	mn = [phs minVal];
	mx -= mn;
	mx = round(mx / 10.0) * 10.0;
//printf("mx  = %f\n", mx);

	return 2 * M_PI / mx;
}

void
mgphs_to_cpx1(NSString *base, NSString *inser, NSString *ser, int nslc)
{
	float		scl = M_PI / 10000;
//	NSString	*inser = @"orig_images";
	NSString	*mgpath, *phspath, *cpxpath;
	RecImage	*img, *phs, *slc, *mimg;
	int			i, nimg, navg;
	RecLoop		*slcLp, *bLp, *avgLp, *yLp, *xLp;

	mgpath = [NSString stringWithFormat:@"%@/%@/%@-mg", base, inser, ser];
	phspath = [NSString stringWithFormat:@"%@/%@/%@-phs", base, inser, ser];

	img = [RecImage imageWithKOImage:mgpath];
	phs = [RecImage imageWithKOImage:phspath];
	nimg = [img zDim];
	navg = nimg / nslc;
	// scale each slice
	for (i = 0; i < nimg; i++) {
		slc = [phs sliceAtIndex:i];
		scl = phs_scale(slc);
		[slc multByConst:scl];
		[phs copySlice:slc atIndex:i];
	}
	[img makeComplexWithPhs:phs];
	[img pcorr];
	[img checkNeg0];

	slcLp = [RecLoop loopWithDataLength:nslc];
	avgLp = [RecLoop loopWithDataLength:navg];
	bLp   = [RecLoop loopWithDataLength:1];
	xLp   = [RecLoop loopWithDataLength:[img xDim]];
	yLp   = [RecLoop loopWithDataLength:[img yDim]];

	mimg = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:avgLp, bLp, slcLp, yLp, xLp, nil];
	[mimg copyImageData:img];

//	cpxpath = [NSString stringWithFormat:@"%@/%@/img.cpx", base, ser];
//	[mimg saveAsKOImage:cpxpath];
	for (i = 0; i < nslc; i++) {
		cpxpath = [NSString stringWithFormat:@"%@/%@/img%d.cpx", base, ser, i];
		slc = [mimg sliceAtIndex:i forLoop:[mimg zLoop]];
		slc = [slc sliceAtIndex:1 forLoop:[slc zLoop]];
		[slc saveAsKOImage:cpxpath];
	}
}

void
mgphs_to_cpx2(NSString *base, NSString *inser, NSString *ser, int nslc)
{
	float		scl = M_PI / 10000;
//	NSString	*inser = @"orig_images";
	NSString	*mgpath, *phspath, *cpxpath;
	RecImage	*img, *phs, *slc, *mimg;
	int			i, nimg, navg;
	RecLoop		*slcLp, *bLp, *avgLp, *yLp, *xLp;

	mgpath = [NSString stringWithFormat:@"%@/%@/%@-mg", base, inser, ser];
	phspath = [NSString stringWithFormat:@"%@/%@/%@-phs", base, inser, ser];

	img = [RecImage imageWithKOImage:mgpath];
	phs = [RecImage imageWithKOImage:phspath];
	nimg = [img zDim];
	navg = nimg / nslc / 2;
	// scale each slice
	for (i = 0; i < nimg; i++) {
		slc = [phs sliceAtIndex:i];
		scl = phs_scale(slc);
		[slc multByConst:scl];
		[phs copySlice:slc atIndex:i];
	}
	[img makeComplexWithPhs:phs];

	// for PC data
	if (1) {
		cpxpath = [NSString stringWithFormat:@"%@/%@/img%d.cpx", base, ser, 0];
		[img pcorr];
		[img saveAsKOImage:cpxpath];
		return;		// for PC data
	}


	[img pcorr];
	[img checkNeg0];

	slcLp = [RecLoop loopWithDataLength:nslc];
	avgLp = [RecLoop loopWithDataLength:navg];
	bLp   = [RecLoop loopWithDataLength:2];
	xLp   = [RecLoop loopWithDataLength:[img xDim]];
	yLp   = [RecLoop loopWithDataLength:[img yDim]];

	mimg = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:avgLp, bLp, slcLp, yLp, xLp, nil];
	[mimg copyImageData:img];

//	cpxpath = [NSString stringWithFormat:@"%@/%@/img.cpx", base, ser];
//	[mimg saveAsKOImage:cpxpath];
	for (i = 0; i < nslc; i++) {
		cpxpath = [NSString stringWithFormat:@"%@/%@/img%d.cpx", base, ser, i];
		slc = [mimg sliceAtIndex:i forLoop:[mimg zLoop]];
		slc = [slc sliceAtIndex:1 forLoop:[slc zLoop]];
		[slc saveAsKOImage:cpxpath];
	}
}

void
reslice()
{
//	NSString	*base = @"../toshiba_images/DWI-nasu-1";
	NSString	*base = @"../toshiba_images/DWI-nasu-2";
//	NSString	*loc = @"Run63070/SSFP_ls_hi";
//	NSString	*loc = @"Run63070/SSFP_ls_lo";
//	NSString	*loc = @"Run63068/SSFP_sag_1";
//	NSString	*loc = @"Run63069/SSFP_sag_1";
	NSString	*loc = @"Run63071/SSFP_sag_1";
	NSString	*path;
	RecImage	*img;

	path = [NSString stringWithFormat:@"%@/%@", base, loc];
	img = [RecImage imageWithKOImage:path];
	[img swapLoop:[img xLoop] withLoop:[img zLoop]];
	path = [NSString stringWithFormat:@"%@/%@.time", base, loc];
	[img saveAsKOImage:path];
}

// read AQP4 file (convert to KOImage)
void
read_raw_aqp()
{
	FILE		*fp;
	char		*path = "/Users/oshio/images/NCI/NIRS/AQP4/AQP4.dat";
	RecImage	*img1, *img2;
	// scan order is: y -> z -> x
	// probably y is centric
	int			xdim = 128;
	int			ydim = 128;
	int			zdim = 14;
	int			i, j, k;
	float		*p1, *q1;
	float		*p2, *q2;
	int			*buf = (int *)malloc(sizeof(int) * xdim * 2);

	img1 = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:xdim yDim:ydim zDim:zdim];
	img2 = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:xdim yDim:ydim zDim:zdim];
	p1 = [img1 real];
	q1 = [img1 imag];
	p2 = [img2 real];
	q2 = [img2 imag];
	fp = fopen(path, "r");
	for (i = 0; i < ydim; i++) {
		for (k = 0; k < zdim; k++) {
			fread(buf, sizeof(int), xdim * 2, fp);
			for (j = 0; j < xdim; j++) {
				p1[k * xdim * ydim + i * xdim + j] = buf[j*2];
				q1[k * xdim * ydim + i * xdim + j] = buf[j*2 + 1];
			}
			fread(buf, sizeof(int), xdim * 2, fp);
			for (j = 0; j < xdim; j++) {
				p2[k * xdim * ydim + i * xdim + j] = buf[j*2];
				q2[k * xdim * ydim + i * xdim + j] = buf[j*2 + 1];
			}
		}
	}
	fclose(fp);
	[img1 fft2d:REC_FORWARD];
	[img1 yFlip];
[img1 crop:[img1 xLoop] to:64 startAt:64];
	[img1 saveAsKOImage:@"/Users/oshio/images/NCI/NIRS/AQP4/AQP4_1.img"];
	[img2 fft2d:REC_FORWARD];
	[img2 yFlip];
img1 = [img2 copy];
[img1 crop:[img1 xLoop] to:64 startAt:0];
	[img1 saveAsKOImage:@"/Users/oshio/images/NCI/NIRS/AQP4/AQP4_2.img"];
img1 = [img2 copy];
[img1 crop:[img1 xLoop] to:64 startAt:64];
	[img1 saveAsKOImage:@"/Users/oshio/images/NCI/NIRS/AQP4/AQP4_3.img"];
}

void
rinkan_1()
{
	RecImage	*img, *avg, *slc;
	NSString	*base = @"/Users/oshio/Projects6/RecKit/toshiba_images/DWI-rinkan-1/Results";
	NSString	*ser = @"SI1";
	NSString	*path;
	int			i, n;
	float		rms;

//	path = [NSString stringWithFormat:@"%@/%@/b200-200-mg", base, ser];
//	path = [NSString stringWithFormat:@"%@/%@/b500-500-mg", base, ser];
	path = [NSString stringWithFormat:@"%@/%@/b1000-1000-mg", base, ser];
	img = [RecImage imageWithKOImage:path];

	n = [img zDim];
	for (i = 0; i < n; i++) {
		slc = [img sliceAtIndex:i];
		rms = [slc rmsVal];
		[slc multByConst:1000/rms];
		[img copySlice:slc atIndex:i];
	}

	avg = [img avgForLoop:[img zLoop]];
	path = [NSString stringWithFormat:@"%@/%@/b200-200-avg", base, ser];
	[avg saveAsKOImage:path];

	[img subImage:avg];
	path = [NSString stringWithFormat:@"%@/%@/b200-200-delta", base, ser];
	[img saveAsKOImage:path];

	mg_pca();

	printf("test\n");

}

void
rinkan_4()
{
	NSString	*base = @"../toshiba_images/DWI-rinkan-4";
	NSString	*inser = @"orig_images";
	
	mgphs_to_cpx2(base,  inser, @"b100-cor-si", 3);
	mgphs_to_cpx2(base,  inser, @"b200-cor-si", 3);
	mgphs_to_cpx2(base,  inser, @"b500-cor-si", 3);
}

void
rinkan_5()
{
	NSString	*base = @"../toshiba_images/DWI-rinkan-5";
	NSString	*inser = @"orig_images";
	
	mgphs_to_cpx2(base,  inser, @"b200-sag-ap-1", 3);
	mgphs_to_cpx2(base,  inser, @"b200-sag-ap-2", 3);
	mgphs_to_cpx2(base,  inser, @"b200-sag-si-1", 3);
	mgphs_to_cpx2(base,  inser, @"b200-sag-si-2", 3);
}

void
nasu_1()
{
	NSString	*base = @"../toshiba_images/DWI-nasu-1/Run62546";
//	NSString	*base = @"../toshiba_images/DWI-nasu-1/Run62547";
//	NSString	*base = @"../toshiba_images/DWI-nasu-1/Run62548";
//	NSString	*base = @"../toshiba_images/DWI-nasu-1/Run62549";
	NSString	*inser = @"tmp";
	

	mgphs_to_cpx1(base, inser, @"b1000-cor-lr", 1);
	mgphs_to_cpx1(base, inser, @"b1000-cor-si", 1);
	mgphs_to_cpx1(base, inser, @"b1000-cor-ap", 1);

//	mgphs_to_cpx1(base, inser, @"b1000-ax-lr", 1);
//	mgphs_to_cpx1(base, inser, @"b1000-ax-si", 1);
//	mgphs_to_cpx1(base, inser, @"b1000-ax-ap", 1);

//	mgphs_to_cpx1(base, inser, @"b200-cor-lr", 1);
//	mgphs_to_cpx1(base, inser, @"b200-cor-si", 1);
//	mgphs_to_cpx1(base, inser, @"b200-cor-ap", 1);

//	mgphs_to_cpx1(base, inser, @"ssfp-2400", 1);
//	mgphs_to_cpx1(base, inser, @"ssfp-2500", 1);
//	mgphs_to_cpx1(base, inser, @"ssfp-2600", 1);
//	mgphs_to_cpx1(base, inser, @"ssfp-2700", 1);
}

void
read_csv()
{
	NSString	*ser = @"../toshiba_images/DWI-nasu-1";
	int			run = 62549;
	int			scan;	// ax, si
	NSString	*path;
	const char	*cPath = [path UTF8String];
	char		*buf;
	int			i, nline, nimg, flg, ofs;
	FILE		*fp;
	RecImage	*img;
	RecLoop		*newX;
	float		*p;
	float		tm, resp, per, ecg1, ecg2, trig, adc;
	float		st, ed;

	buf = (char *)malloc(256);
	for (scan = 7000; scan < 30000; scan += 1000) {
		path = [NSString stringWithFormat:@"%@/DVD20180511exp/TEST%d/Run%d.-5004_%d.csv", ser, run - 62544, run, scan];
		cPath = [path UTF8String];
		fp = fopen(cPath, "r");
		if (fp == NULL) continue;


		for (i = 0; ; i++) {
			if (fgets(buf, 256, fp) == NULL) break;
		}
		nline = i;
		printf("%d lines read\n", i);

		img = [RecImage imageOfType:RECIMAGE_REAL xDim:nline yDim:4];
		p = [img data];
		rewind(fp);
		for (i = 0; i < nline; i++) {
			if (fgets(buf, 256, fp) == NULL) break;
			if (i == 0) continue;	// skip format line
			sscanf(buf, "%f,%f,%f,%f,%f,%f,%f", &tm, &resp, &per, &ecg1, &ecg2, &trig, &adc);
		//	printf("%6.3f %6.3f %6.3f %6.3f\n", tm, resp, per, adc * 8000);
			p[i] = tm;
			p[nline * 1 + i] = resp;
			p[nline * 2 + i] = per;
			p[nline * 3 + i] = adc * 5000;
		}
		
		// count adc trig
		nimg = flg = ofs = 0;
		st = ed = 0;
		for (i = 0; i < nline; i++) {
			tm = p[i];
			adc = p[nline * 3 + i];
			if (flg == 0 && adc) {
				nimg++;
				flg = 1;
				if (st == 0) {
					st = tm;
					ofs = i;
				}
			}
			if (flg == 1 && adc == 0) {
				flg = 0;
				ed = tm;
			}
		}
		printf("nimg = %d, tm = %f\n", nimg, ed - st);
	//printf("ofs = %d\n", ofs);
		nimg = 0;
		for (i = 0; i < nline; i++) {
			tm = p[i];
			if (tm < st) continue;
			nimg++;
			if (tm >= ed) {
				break;
			}
		}
	//	printf("nimg = %d, st/ed = %f/%f\n", nimg, st, ed);

	// crop
		newX = [RecLoop loopWithDataLength:nimg];
		[img replaceLoop:[img xLoop] withLoop:newX offset:ofs];
	// sub-sample
		img = [img scale1dLoop:[img xLoop] to:256];

		path = [NSString stringWithFormat:@"%@/Run%d/tmp/phys-%d.img", ser, run, scan];
		[img saveAsKOImage:path];

		fclose(fp);
	}
	free(buf);
}

void
read_pc()
{
	NSString			*base = @"../toshiba_images/DWI-nasu-5";
	NSString			*inser = @"OrigImages";

	mgphs_to_cpx2(base, inser, @"3V", 1);
	mgphs_to_cpx2(base, inser, @"3Vu", 1);
	mgphs_to_cpx2(base, inser, @"4V", 1);
	mgphs_to_cpx2(base, inser, @"4Vu", 1);
}

void
pc_roi()
{
	NSString			*base = @"../toshiba_images/DWI-nasu-5";
	NSString			*inser = @"4V";
	NSString			*path;
	int					cx = 63, cy = 75;
	float				rr = 5;
	float				r, x, y;
	int					rx, ry;
	int					siz = 12;
	int					i, j, ix, dim;
	int					k, nImg;
	RecLoopControl		*lc;
	RecImage			*img;
	float				*p;
	float				outer, inner, accum;
	int					no, ni;
	BOOL				rect, circ;

	path = 	[NSString stringWithFormat:@"%@/%@/img0.cpx", base, inser];
	img = [RecImage imageWithKOImage:path];
	[img phase];
	dim = [img xDim];	// square
	nImg = [img zDim];

	lc = [img control];
	[lc deactivateXY];
	rx = cx - siz/2;
	ry = cy - siz/2;
	accum = 0;
	for (k = 0; k < nImg; k++) {
		p = [img currentDataWithControl:lc];
		// first pass
		no = ni = 0;
		outer = inner = 0.0;
		for (i = 0; i < dim; i++) {
			for (j = 0; j < dim; j++) {
				ix = i * dim + j;
				r = sqrt((i - cy) * (i - cy) + (j - cx) * (j - cx));
				rect = circ = NO;
				if (r < rr) circ = YES;
				if (j > rx && j < rx + siz &&
					i > ry && i < ry + siz) {
					rect = YES;
				}
				if (rect && ! circ) {
					outer += p[ix];
					no += 1;
				//	p[ix] += 0.2;
				}
				if (circ) {
				//	p[ix] -= 0.2;
				}
			}
		}
		outer /= no;

		// second pass
		ni = no = 0;
		inner = 0;
		for (i = 0; i < dim; i++) {
			for (j = 0; j < dim; j++) {
				ix = i * dim + j;
				r = sqrt((i - cy) * (i - cy) + (j - cx) * (j - cx));
				rect = circ = NO;
				if (r < rr) circ = YES;
				if (j > rx && j < rx + siz &&
					i > ry && i < ry + siz) {
					rect = YES;
				}
				if (rect) {
					p[ix] -= outer;
				//	p[ix] += 0.2;
				}
			//	if (circ) {
				if (rect) {
					inner += p[ix];
					ni += 1;
				//	p[ix] -= 0.2;
				}
			}
		}
		inner /= ni;
		accum += inner;
		
		printf("%d %f %f\n", k, inner, accum);

		[lc increment];
	}

	path = 	[NSString stringWithFormat:@"%@/%@/phase.roi", base, inser];
	[img saveAsKOImage:path];
}

void
pc_avg()
{
	NSString	*base = @"../toshiba_images/DWI-nasu-5";
//	NSString	*inser = @"3V";
//	NSString	*inser = @"3Vu";
//	NSString	*inser = @"4V";
	NSString	*inser = @"4Vu";
	NSString	*path;
	RecImage	*img, *avg;
	int			dim, nImg;

	path = 	[NSString stringWithFormat:@"%@/%@/img0.cpx", base, inser];
	img = [RecImage imageWithKOImage:path];
	[img phase];
	[img removeSliceAtIndex:0 forLoop:[img zLoop]];
	path = 	[NSString stringWithFormat:@"%@/%@/img0.phs", base, inser];
	[img saveAsKOImage:path];

	avg = [img avgForLoop:[img zLoop]];
	[avg gauss2DHP:0.05 frac:1.0];
	path = 	[NSString stringWithFormat:@"%@/%@/phase.avg", base, inser];
	[avg saveAsKOImage:path];

	avg = [img copy];

	[avg gauss1DLP:0.1 forLoop:[avg zLoop]];
    [avg gauss2DHP:0.05 frac:1.0];
	path = 	[NSString stringWithFormat:@"%@/%@/phase.mvavg", base, inser];
	[avg saveAsKOImage:path];
	
}

// mg PCA
void
mg_pca()
{
//	NSString			*path = @"../toshiba_images/DWI-rinkan-4/b200-cor-si/img1.cpx";
//	NSString			*path = @"../toshiba_images/DWI-rinkan-2/b1000-proc/b1000.img";
//========================
	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62546/b1000-cor-si";
//	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62546/b1000-cor-ap";
//	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62546/b1000-cor-lr";
	NSString			*path;
	RecImage			*img;
	RecImage			*mask;
	float				*p;
	int					i, j, xDim, yDim;

	Num_svd_result		*sres;
//	Num_ica_result		*ires;
	Num_mat				*A;
	RecImage			*a_im, *u_im, *vt_im;
	RecImage			*avg, *frq;
	int					np, n;

	path = 	[NSString stringWithFormat:@"%@/img0.cpx", base];
	img = [RecImage imageWithKOImage:path];
	avg = [img avgForLoop:[img zLoop]];
	path = 	[NSString stringWithFormat:@"%@/IMG_avg", base];
	[avg saveAsKOImage:path];
	[img magnitude];

	avg = [img varForLoop:[img zLoop]];
	path = 	[NSString stringWithFormat:@"%@/IMG_var", base];
	[avg saveAsKOImage:path];

	avg = [img avgForLoop:[img zLoop]];
	[img subImage:avg];

// === freq analysis
	frq = [img copy];
	[frq fft1d:[frq zLoop] direction:REC_INVERSE];
	[frq gauss1DLP:0.2 forLoop:[frq zLoop]];
	[frq gauss2DLP:0.4];
	path = 	[NSString stringWithFormat:@"%@/IMG_mg_f", base];
	[frq saveAsKOImage:path];

//exit(0);

// === make mask
	xDim = [img xDim];
	yDim = [img yDim];
	mask = [RecImage imageOfType:RECIMAGE_REAL withLoops:[img yLoop], [img xLoop], nil];
	p = [mask data];
	for (i = 0; i < yDim; i++) {
		for (j = 0; j < xDim; j++) {
			if (i > yDim * 0.25 && i < yDim * 0.75 &&
				j > xDim * 0.3 && j < xDim * 0.7) {
				p[i * xDim + j] = 1;
			}
		}
	}
	[img multByImage:mask];
	path = [NSString stringWithFormat:@"%@/IMG_mask", base];
	[mask saveAsKOImage:path];
	path = [NSString stringWithFormat:@"%@/IMG_img_m", base];
	[img saveAsKOImage:path];

	np = [img xDim] * [img yDim];
	n = [img zDim];
	a_im = [RecImage imageOfType:RECIMAGE_REAL xDim:np yDim:n];
//	a_im = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:np yDim:n];
	[a_im copyImageData:img];
	A = Num_im_to_m(a_im);
	sres = Num_svd(A);
	u_im = Num_m_to_im(sres->U);
	[u_im trans];
//	Num_scale_rows(sres->Vt, sres->s);
	vt_im = Num_m_to_im(sres->Vt);

	[u_im crop:[u_im yLoop] to:20 startAt:0];
	path = [NSString stringWithFormat:@"%@/IMG_U", base];
	[u_im saveAsKOImage:path];
//	[img copyImageData:vt_im];

	path = [NSString stringWithFormat:@"%@/IMG_PCA", base];
	[img saveAsKOImage:path];
	printf("end\n");
}

// phs PCA
void
phs_pca()
{
//====
//	NSString			*path = @"../toshiba_images/DWI-rinkan-4/b200-cor-si/img1.cpx";
//	NSString			*path = @"../toshiba_images/DWI-rinkan-5/b200-sag-si-1/img1.cpx";
//	NSString			*path = @"../toshiba_images/DWI-rinkan-2/b1000-proc/b1000.img";
//========================
	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62546/b1000-cor-si";
//	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62546/b1000-cor-ap";
//	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62546/b1000-cor-lr";
	NSString			*path;
	RecImage			*img;

	Num_svd_result		*sres;
	Num_ica_result		*ires;
	Num_mat				*A;
	RecImage			*a_im, *u_im, *vt_im;
	RecImage			*avg, *phs, *mag, *frq;
	RecImage			*mask;
	RecImage			*eimg, *tmp;
	int					i, np, n;
	int					orig_dim;
	int					ncomp = 10;

	path = 	[NSString stringWithFormat:@"%@/img0.cpx", base];

	img = [RecImage imageWithKOImage:path];
	path = 	[NSString stringWithFormat:@"%@/IMG_in", base];
	[img saveAsKOImage:path];

// pre-proc
	avg = [img avgForLoop:[img zLoop]];
	path = 	[NSString stringWithFormat:@"%@/IMG_avg", base];
	[avg saveAsKOImage:path];
	mask = [avg copy];
//	[mask thresAt:1.0e-8];
	[mask thresAt:0.05];
	path = 	[NSString stringWithFormat:@"%@/IMG_mask", base];
	[mask saveAsKOImage:path];

	[img multByImage:mask];
	[img phase];
	phs = [img unwrapEst2d];

	[phs multByImage:mask];
	avg = [phs avgForLoop:[phs zLoop]];
	[phs subImage:avg];
	path = 	[NSString stringWithFormat:@"%@/IMG_unwrap", base];
	[phs saveAsKOImage:path];

// === freq analysis
	frq = [phs copy];
	[frq fft1d:[frq zLoop] direction:REC_INVERSE];
	[frq gauss1DLP:0.2 forLoop:[frq zLoop]];
	path = 	[NSString stringWithFormat:@"%@/IMG_phs_f", base];
	[frq saveAsKOImage:path];

exit(0);

// === PCA -> make func
	np = [phs xDim] * [phs yDim];
	n = [phs zDim];
	a_im = [RecImage imageOfType:RECIMAGE_REAL xDim:np yDim:n];
//	a_im = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:np yDim:n];
	[a_im copyImageData:phs];
	A = Num_im_to_m(a_im);
	sres = Num_svd(A);
	u_im = Num_m_to_im(sres->U);
	[u_im trans];
	vt_im = Num_m_to_im(sres->Vt);
	path = @"U_img.img";
	[u_im saveAsKOImage:path];
	[phs copyImageData:vt_im];
	path = @"E_phs.img";
	[phs saveAsKOImage:path];
	for (i = 0; i < sres->s->n; i++) {
		printf("%d %f\n", i, sres->s->data[i]);
	}

// === mult mag image with U
	if (0) {
		[a_im copyImageData:mag];
		A =  Num_im_to_m(a_im);
		Num_clear_mat(sres->Vt);
		Num_mtmul(sres->Vt, sres->U, YES, A, NO);
		vt_im = Num_m_to_im(sres->Vt);
		[mag copyImageData:vt_im];
		[mag multByImage:mask];
		[mag saveAsKOImage:@"E_mag.img"];
	}

exit(0);

// === ICA
	ires = Num_ica(A, ncomp);
	tmp = Num_m_to_im(ires->WK);
	eimg = [RecImage imageOfType:RECIMAGE_REAL xDim:[img xDim] yDim:[img yDim] zDim:ncomp];
	[eimg copyImageData:tmp];
	[eimg saveAsKOImage:@"IMG_map"];
	saveAsKOImage(ires->WX, @"IMG_out");

	Num_free_svd_result(sres);
	Num_free_ica_result(ires);

}


// cpx PCA
void
cpx_pca()
{
//	NSString			*path = @"../toshiba_images/DWI-rinkan-4/b200-cor-si/img1.cpx";
//	NSString			*path = @"../toshiba_images/DWI-rinkan-5/b200-sag-si-1/img1.cpx";
	NSString			*path = @"../toshiba_images/DWI-rinkan-2/b1000-proc/b1000.img";
	RecImage			*img;

	Num_svd_result		*sres;
	Num_ica_result		*ires;
	Num_mat				*A;
	RecImage			*a_im, *u_im, *vt_im;
	RecImage			*avg, *mg, *phs;
	RecImage			*mask;
	RecImage			*eimg, *tmp;
	int					i, j, np, n, ix;
	float				*p, *q, *pp;
	int					orig_dim;
	int					ncomp = 10;

	img = [RecImage imageWithKOImage:path];
	orig_dim = [img xDim];
//	[img zeroFillToPo2];
	[img cropToPo2];
[img saveAsKOImage:@"IMG_in"];

// mask
	avg = [img avgForLoop:[img zLoop]];
	[avg saveAsKOImage:@"IMG_avg"];
	mask = [avg copy];
//	[mask thresAt:1.0e-8];
	[mask thresAt:0.05];
	[mask saveAsKOImage:@"IMG_mask"];

	[img multByImage:mask];

	mg = [img copy];
	[mg magnitude];
	avg = [mg avgForLoop:[mg zLoop]];
	[mg subImage:avg];
[mg saveAsKOImage:@"IMG_mg"];

//	phs = [img unwrap2d];
	phs = [img copy];
	[phs phase];
	phs = [phs unwrapEst2d];
//	[img crop:[img xLoop] to:orig_dim];
//	[img crop:[img yLoop] to:orig_dim];

	[phs multByImage:mask];
	avg = [phs avgForLoop:[phs zLoop]];
	[phs subImage:avg];
	[phs multByConst:100];
	[phs saveAsKOImage:@"IMG_phs"];

	img = [mg copy];
	[img makeComplexWithIm:phs];
	[img saveAsKOImage:@"IMG_cpx"];

// === PCA -> make func
	np = [img xDim] * [phs yDim] * 2;
	n = [img zDim];
	a_im = [RecImage imageOfType:RECIMAGE_REAL xDim:np yDim:n];
//	[a_im copyImageData:img];
//	[a_im saveAsKOImage:@"IMG_A"];
//exit(0);

	p = [img real];	// dst
	q = [img imag];
	pp = [a_im data];		// src

	ix = 0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < np/2; j++) {
			pp[ix++] = p[i * np/2 + j];
		}
		for (j = 0; j < np/2; j++) {
			pp[ix++] = q[i * np/2 + j];
		}
	}
printf("%d %d\n", ix, [a_im dataLength]);
[a_im dumpLoops];

//	[a_im saveAsKOImage:@"IMG_A"];
//exit(0);



	A = Num_im_to_m(a_im);
	sres = Num_svd(A);
	u_im = Num_m_to_im(sres->U);
	[u_im trans];
	vt_im = Num_m_to_im(sres->Vt);
	path = @"U_img.img";
	[u_im saveAsKOImage:path];
//	[img copyImageData:vt_im];
	pp = [vt_im data];
	p = [img real];
	q = [img imag];
	ix = 0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < np/2; j++) {
			p[i * np/2 + j] = pp[ix++];
		}
		for (j = 0; j < np/2; j++) {
			q[i * np/2 + j] = pp[ix++];
		}
	}


	path = @"E_img.img";
	[img saveAsKOImage:path];
	for (i = 0; i < sres->s->n; i++) {
		printf("%d %f\n", i, sres->s->data[i]);
	}

// === ICA
	ires = Num_ica(A, ncomp);
	tmp = Num_m_to_im(ires->WK);
	eimg = [RecImage imageOfType:RECIMAGE_REAL xDim:[img xDim] yDim:[img yDim] zDim:ncomp];
	[eimg copyImageData:tmp];
	[eimg saveAsKOImage:@"IMG_map"];
	saveAsKOImage(ires->WX, @"IMG_out");

	Num_free_svd_result(sres);
	Num_free_ica_result(ires);

}

void
mp_freq()
{
// ### make these NSArray
//	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62546/b1000-ax-si";
//	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62546/b1000-ax-lr";
//	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62546/b1000-ax-ap";

//	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62546/b1000-cor-si";
//	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62546/b1000-cor-lr";
	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62546/b1000-cor-ap";

//	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62547/b1000-cor-si";
//	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62547/b1000-cor-lr";
//	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62547/b1000-cor-ap";

//	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62547/b200-cor-si";
//	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62547/b200-cor-lr";
//	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62547/b200-cor-ap";

//	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62548/b1000-cor-si";
//	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62548/b1000-cor-lr";
//	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62548/b1000-cor-ap";

//	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62548/b200-cor-si";
//	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62548/b200-cor-lr";
//	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62548/b200-cor-ap";

//	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62549/b1000-cor-si";
//	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62549/b1000-cor-lr";
//	NSString			*base = @"../toshiba_images/DWI-nasu-1/Run62549/b1000-cor-ap";

	NSString			*path;
	RecImage			*img;

	RecImage			*avg, *phs, *mag, *frq;
	RecImage			*mask;

	Num_svd_result		*sres;
//	Num_ica_result		*ires;
	Num_mat				*A;
	RecImage			*a_im, *u_im, *vt_im;
	int					i, np, n;

	path = 	[NSString stringWithFormat:@"%@/img0.cpx", base];
	img = [RecImage imageWithKOImage:path];

// pre-proc
	avg = [img avgForLoop:[img zLoop]];
	path = 	[NSString stringWithFormat:@"%@/IMG_avg", base];
	[avg saveAsKOImage:path];
	mask = [avg copy];
//	[mask thresAt:1.0e-8];
	[mask thresAt:0.05];
	path = 	[NSString stringWithFormat:@"%@/IMG_mask", base];
	[mask saveAsKOImage:path];

// mag
	mag = [img copy];
	[mag magnitude];
	avg = [mag avgForLoop:[img zLoop]];
	[mag subImage:avg];
	path = 	[NSString stringWithFormat:@"%@/IMG_mag_0", base];
	[mag saveAsKOImage:path];

	frq = [mag copy];
	[frq fft1d:[img zLoop] direction:REC_INVERSE];
	path = 	[NSString stringWithFormat:@"%@/IMG_mg_f", base];
	[frq saveAsKOImage:path];

// phase
	[img multByImage:mask];
	[img phase];
	phs = [img unwrapEst2d];
	[phs multByImage:mask];
	avg = [phs avgForLoop:[phs zLoop]];
	[phs subImage:avg];
	path = 	[NSString stringWithFormat:@"%@/IMG_unwrap", base];
	[phs saveAsKOImage:path];

// === freq analysis
	frq = [phs copy];
	[frq fft1d:[frq zLoop] direction:REC_INVERSE];
//	[frq gauss1DLP:0.2 forLoop:[frq zLoop]];
	path = 	[NSString stringWithFormat:@"%@/IMG_phs_f", base];
	[frq saveAsKOImage:path];

// ==== PCA phs (## add multi-band later)
	np = [phs xDim] * [phs yDim];
	n = [phs zDim];
	a_im = [RecImage imageOfType:RECIMAGE_REAL xDim:np yDim:n];
	[a_im copyImageData:phs];
	A = Num_im_to_m(a_im);
	sres = Num_svd(A);
	u_im = Num_m_to_im(sres->U);
	[u_im trans];
	[u_im crop:[u_im yLoop] to:10 startAt:0];
	vt_im = Num_m_to_im(sres->Vt);
	path = 	[NSString stringWithFormat:@"%@/IMG_phs_U", base];
	[u_im saveAsKOImage:path];
	[phs copyImageData:vt_im];
	path = 	[NSString stringWithFormat:@"%@/IMG_phs_E", base];
	[phs saveAsKOImage:path];
//	for (i = 0; i < sres->s->n; i++) {
//		printf("%d %f\n", i, sres->s->data[i]);
//	}

// ==== PCA mg (## add multi-band later)
	n = [mag zDim] - 10;
	np = [mag xDim] * [mag yDim];
	a_im = [RecImage imageOfType:RECIMAGE_REAL xDim:np yDim:n];
	[mag crop:[mag zLoop] to:n startAt:10];
	path = 	[NSString stringWithFormat:@"%@/IMG_mag_1", base];
	[mag saveAsKOImage:path];
	[a_im copyImageData:mag];
	A = Num_im_to_m(a_im);
	sres = Num_svd(A);
	u_im = Num_m_to_im(sres->U);
	[u_im trans];
	[u_im crop:[u_im yLoop] to:10 startAt:0];
	vt_im = Num_m_to_im(sres->Vt);
	path = 	[NSString stringWithFormat:@"%@/IMG_mg_U", base];
	[u_im saveAsKOImage:path];
	[mag copyImageData:vt_im];
	path = 	[NSString stringWithFormat:@"%@/IMG_mg_E", base];
	[mag saveAsKOImage:path];
//	for (i = 0; i < sres->s->n; i++) {
//		printf("%d %f\n", i, sres->s->data[i]);
//	}


	
	Num_free_svd_result(sres);
}

// cardiac phase avg
void
card_avg()
{
	NSString		*base = @"../toshiba_images/DWI-nasu-1/Run62546/b1000-cor-si";
	NSString		*path;
	RecImage		*phys, *unwrap, *tab, *pulse, *avg, *ref, *corr;
	float			sft[256];
	RecImage		*img_p;			// image pulase data (1 x tDim2)
	RecImage		*slc, *corSlc;
	int				i, j, pk;
	int				nPhs = 8;
//	int				cnt[50];
	int				tDim1, tDim2, nPulse;
	int				xDim, yDim;
	int				ix1, ix2;
	float			*p, mn, mx, val, del;
	float			*q;
	float			*tabt, *tabd;

	printf("card avg\n");

	// read phys, unwrap
	path   = [NSString stringWithFormat:@"%@/phys.img", base];
	phys   = [RecImage imageWithKOImage:path];
	path   = [NSString stringWithFormat:@"%@/IMG_unwrap", base];
	unwrap = [RecImage imageWithKOImage:path];
	
	// make table img
	tDim1 = [phys xDim];
	xDim = [unwrap xDim];
	yDim = [unwrap yDim];
	tDim2 = [unwrap zDim];
	tab = [RecImage imageOfType:RECIMAGE_REAL xDim:256 yDim:2];	// peak pos, duration
	tabt = [tab data];
	tabd = tabt + [tab dataLength];

	phys = [phys sliceAtIndex:2 forLoop:[phys yLoop]];	// 0: time, 1: resp, 2: card

	// copy image phase array (1D)
	img_p = [RecImage imageOfType:RECIMAGE_REAL xDim:tDim2];
	p = [unwrap data];
	q = [img_p data];
	for (i = 0; i < tDim2; i++) {
		q[i] = p[i * xDim * yDim + yDim * 31 + 31];
	//	printf("%d %f\n", i, q[i]);
	}

	// make tab
	p = [phys data];
	mn = [phys meanVal];
	mx = [phys maxVal];
	mn = (mn*2 + mx) / 3;
	[phys addConst:-mn];

	pk = 0;
	nPulse = 0;
	for (i = 0; i < tDim1; i++) {
		val = p[i];
		if (val > 0) {
			for (j = i; j < tDim1; j++) {
				if (p[j] < val) break;
				val = p[j];
			}
		//	printf("%d %8.1f %d\n", j, p[j], j - pk);
			tabt[nPulse] = j * 0.01;	// sec
			tabd[nPulse++] = (j - pk) * 0.01;
			pk = j;
			for (; i < tDim1; i++) {
				if (p[i] < 0) break;
			}
		}
	}
	del = 0; //-0.7;
	for (i = 0; i < nPulse; i++) {
		tabt[i] += del;
	//	printf("%d %8.2f %8.2f\n", i, tabt[i], tabd[i]);
	}

printf("nPulse = %d\n", nPulse);

//===== break into  pulses
	pulse = [RecImage imageOfType:RECIMAGE_REAL xDim:nPhs yDim:nPulse];
	q = [pulse data];
	p = [img_p data];
	for (i = 0; i < nPulse; i++) {
		ix1 = (int)(nPhs * i);
		ix2 = (int)((tabt[i]) * 5.05) - 5;
		for (j = 0; j < nPhs; j++) {
			q[ix1] = p[ix2];	// ###
			ix1 += 1;
			ix2 += 1;
		}
	}
	path	= [NSString stringWithFormat:@"%@/IMG_pulse.img", base];
	[pulse saveAsKOImage:path];

// ===== iterative correlation / averaging proc
	avg = [pulse copy];
	avg = [avg scale1dLoop:[avg xLoop] to:[avg xDim] * 2];
	path	= [NSString stringWithFormat:@"%@/IMG_pavg.img", base];
	[avg saveAsKOImage:path];
	ref = [avg avgForLoop:[avg yLoop]];
	path	= [NSString stringWithFormat:@"%@/IMG_pref.img", base];
	[ref saveAsKOImage:path];

	corr = [RecImage imageWithImage:avg];
	// iteration step
	for (i = 0; i < [avg yDim]; i++) {
		float *p;
		slc = [avg sliceAtIndex:i forLoop:[avg yLoop]];
		corSlc = [slc xCorrelationWith:ref width:0.2 triFilt:YES];
		p = [corSlc data];
		sft[i] = Rec_find_peak(p, 1, [corSlc xDim]);
	//	printf("%d %f\n", i, sft[i] - [corSlc xDim]/2);
		[corr copySlice:corSlc atIndex:i forLoop:[avg yLoop]];
	}
	path	= [NSString stringWithFormat:@"%@/IMG_corr.img", base];
	[corr saveAsKOImage:path];
	
	// avg with shift
	p = [pulse data];
	for (i = 0; i < [pulse yDim]; i++) {
		for (j = 0; j < [pulse xDim]; j++) {
			printf("%f %5.1f\n", (float)j - (sft[i] - [pulse xDim]) * 0.5 + 0.5, p[i * [pulse xDim] + j]);
		//	printf("%f %5.1f\n", (float)j + (tabt[i] - (int)tabt[i]), p[i * [pulse xDim] + j]);
		}
		printf("!eoc\n");	// not working
	}
exit(0);
		
	
}

void
card_avg2()
{
}

void
t2sim()
{
	int		i;
	float	t2br = 110, dbr = 1.0;				// brain (gray)
	float	t2b = 100,  db  = 5.0, fb = 0.05;	// blood
	float	t2c = 2000, dc  = 3.0, dci = 7.0;	// CSF / ISF
	float	fi = 0.15;
	float	sbr, sb, sc, si;
	float	t;
	float	b = 0.2;

	printf("T2 sim\n");

	b = 0.0;
	for (i = 0; i <= 40; i++) {
		t = 30.0 * i;
		sbr = 1.0 * exp(-b * dbr) * exp(-t / t2br);		// brain (gray)
		sb  = fb  * exp(-b * db)  * exp(-t / t2b);		// blood
		si  = fi  * exp(-b * dci)  * exp(-t / t2c);		// ISF
		sc  = 1.0  * exp(-b * dc)  * exp(-t / t2c);		// CSF
		printf("%f %f %f %f %f\n", t, sbr, sb, sc, si);
	}

	printf("! ===========\n\n");

	printf("IVIM\n");
	t = 200;
	for (i = 0; i <= 10; i++) {
		b = 0.1 * i;
		sbr = 1.0 * exp(-b * dbr) * exp(-t / t2br);		// blood
		sb  = fb  * exp(-b * db)  * exp(-t / t2b);		// brain (gray)
		si  = fi  * exp(-b * dci)  * exp(-t / t2c);		// ISF
		sc  = 1.0  * exp(-b * dc)  * exp(-t / t2c);	// CSF
		printf("%f %f %f %f %f %f\n", b * 1000, sbr, sb, si, sc, sbr + si);
	}

}

// -> TE = 200 or 500, b = 0 - 400

void
calc_div()
{
	NSString	*base = @"../toshiba_images/DWI-nasu-1";
//	NSString	*loc = @"Run62546/b1000-cor-si";
	NSString	*loc = @"Run62546/b1000-cor-lr";
	BOOL		si = NO;
	NSString	*path;
	RecImage	*img;
	RecImage	*div;
	RecImage	*mask;
	float		*q, *p;
	int			i, j, k, xDim, yDim, zDim;
	

	printf("div test\n");

	path = [NSString stringWithFormat:@"%@/%@/img0.cpx", base, loc];
	img = [RecImage imageWithKOImage:path];
	mask = [img copy];
	mask = [mask avgForLoop:[mask zLoop]];
	[mask magnitude];
	[mask thresAt:0.05];
	path = 	[NSString stringWithFormat:@"%@/%@/IMG_mask", base, loc];
	[mask saveAsKOImage:path];

	[img phase];
	p = [img data];
	div = [img copy];
	q = [div data];
	xDim = [img xDim];
	yDim = [img yDim];
	zDim = [img zDim];

	for (k = 0; k < zDim; k++) {
		if (si) {
			for (i = 1; i < yDim; i++) {
				for (j = 0; j < xDim; j++) {
						q[(k * yDim + i) * xDim + j] = 
							p[(k * yDim + i) * xDim + j] - p[(k * yDim + i - 1) * xDim + j];
				}
			}
		} else {
			for (i = 0; i < yDim; i++) {
				for (j = 1; j < xDim; j++) {
						q[(k * yDim + i) * xDim + j] = 
							p[(k * yDim + i) * xDim + j] - p[(k * yDim + i - 1) * xDim + j-1];
				}
			}
		}
	}
	[div unwrap0d];
	[div multByImage:mask];

	path = [NSString stringWithFormat:@"%@/%@/IMG_div_si", base, loc];
	[div saveAsKOImage:path];
}

void
ivcm()
{
	NSString	*base = @"../toshiba_images/DWI-nasu-9";
	NSString	*path;
	RecImage	*raw, *img, *img_p, *img_m;
	float		*p, *q;
	int			i, j, k, xDim, yDim, zDim, nAx;
	RecLoop		*slc, *sign, *axis, *xLp, *yLp, *avgLp;
	

	printf("IVCM\n");

// ser 1
	path = [NSString stringWithFormat:@"%@/DWI1/img0.mag", base];
	raw = [RecImage imageWithKOImage:path];
	[raw dumpLoops];

	avgLp = [RecLoop loopWithDataLength:4];
	sign = [RecLoop loopWithDataLength:2];	// +, -
	axis = [RecLoop loopWithDataLength:3];	// ss, ro, pe
	slc = [RecLoop loopWithDataLength:42];
	xLp = [RecLoop loopWithDataLength:160];
	yLp = [RecLoop loopWithDataLength:160];
	img = [RecImage imageOfType:RECIMAGE_REAL withLoops:avgLp, axis, sign, slc, yLp, xLp, nil];

	[img copyImageData:raw];
	img = [img avgForLoop:avgLp];
	[img dumpLoops];

	path = [NSString stringWithFormat:@"%@/DWI1/IMG.avg", base];
	[img saveAsKOImage:path];
	
	img_p = [img sliceAtIndex:0 forLoop:sign];
	path = [NSString stringWithFormat:@"%@/DWI1/IMG.pls", base];
//	[img_p saveAsKOImage:path];
	img_m = [img sliceAtIndex:1 forLoop:sign];
	path = [NSString stringWithFormat:@"%@/DWI1/IMG.mns", base];
//	[img_m saveAsKOImage:path];

	// sub
	[img_p subImage:img_m];
	path = [NSString stringWithFormat:@"%@/DWI1/IMG.sub", base];
	[img_p saveAsKOImage:path];

// ser 2
	path = [NSString stringWithFormat:@"%@/DWI2/img0.mag", base];
	raw = [RecImage imageWithKOImage:path];
	[raw dumpLoops];

	avgLp = [RecLoop loopWithDataLength:12];
	sign = [RecLoop loopWithDataLength:2];	// +, -
	axis = [RecLoop loopWithDataLength:1];	// ss, ro, pe
	slc = [RecLoop loopWithDataLength:34];
	xLp = [RecLoop loopWithDataLength:160];
	yLp = [RecLoop loopWithDataLength:160];
	img = [RecImage imageOfType:RECIMAGE_REAL withLoops:avgLp, axis, sign, slc, yLp, xLp, nil];
	[img copyImageData:raw];

[img swapLoop:avgLp withLoop:slc];	// slc, axis, sign, avg
[img swapLoop:avgLp withLoop:sign];	// slc, axis, avg, sign
img = [img avgForLoop:avgLp];
path = [NSString stringWithFormat:@"%@/DWI2/IMG.reorder", base];
	[img saveAsKOImage:path];


	img = [img avgForLoop:avgLp];
	[img dumpLoops];

	path = [NSString stringWithFormat:@"%@/DWI2/IMG.avg", base];
	[img saveAsKOImage:path];
	img_p = [img sliceAtIndex:0 forLoop:sign];
	img_m = [img sliceAtIndex:1 forLoop:sign];

	// sub
	[img_p subImage:img_m];
	path = [NSString stringWithFormat:@"%@/DWI2/IMG.sub", base];
	[img_p saveAsKOImage:path];
}

// quantitative time-SLIP
// started on 5-16-2020
void    ts_1()
{
    NSString    *base = @"../toshiba_images/TS-nasu-1/5s";
//    NSString    *base = @"../toshiba_images/TS-nasu-1/2s";
    NSString    *path;
    RecImage    *img1, *img2, *ref1, *ref2;
    RecImage    *sub1, *sub2;
    RecImage    *roi1, *roi2;
    RecImage    *tmp;
    int         i, j, n, tLen;
    float       *p, sum1, sum2;
    float       *ac1, *ac2;
    BOOL        rot = YES;
    BOOL        removeRef = YES;

    system("rm IMG_*.*");

    path = [NSString stringWithFormat:@"%@/IMG_img1.img", base];
    img1 = [RecImage imageWithKOImage:path];
    path = [NSString stringWithFormat:@"%@/IMG_img2.img", base];
    img2 = [RecImage imageWithKOImage:path];
    path = [NSString stringWithFormat:@"%@/IMG_ref1.img", base];
    ref1 = [RecImage imageWithKOImage:path];
    path = [NSString stringWithFormat:@"%@/IMG_ref2.img", base];
    ref2 = [RecImage imageWithKOImage:path];
    [img2 copyLoopsOf:img1];
    [ref1 copyLoopsOf:img1];
    [ref2 copyLoopsOf:img1];

    if (rot) {
        img1 = [img1 rotByTheta:30.0 * M_PI / 180.0];
        img2 = [img2 rotByTheta:30.0 * M_PI / 180.0];
        ref1 = [ref1 rotByTheta:30.0 * M_PI / 180.0];
        ref2 = [ref2 rotByTheta:30.0 * M_PI / 180.0];
    }
    tLen = [img1 zDim]/2 - 1;
    if (removeRef) {
        [img1 crop:[img1 zLoop] to:tLen startAt:1];
        [img2 crop:[img2 zLoop] to:tLen startAt:1];
    }
    ref1 = [ref1 avgForLoop:[ref1 zLoop]];
    ref2 = [ref2 avgForLoop:[ref2 zLoop]];

// checking subtraction images ... original images are not changed
// 4->3, pos
    tmp = [img1 copy];
//    [ref1 multByConst:0.9];
    [tmp subImage:ref1];
    [tmp saveAsKOImage:@"IMG_3p.img"];
// 3-4, neg
    tmp = [img1 copy];
    [tmp subImage:ref2];
    [tmp negate];
    [tmp saveAsKOImage:@"IMG_4n.img"];

// 3-4, pos
    tmp = [img2 copy];
    [tmp subImage:ref1];
    [tmp saveAsKOImage:@"IMG_4p.img"];

// 4-3, neg
    tmp = [img2 copy];
    [tmp subImage:ref2];
    [tmp negate];
    [tmp saveAsKOImage:@"IMG_3n.img"];

// ROI
    ac1 = (float *)malloc(sizeof(float) * tLen);
    ac2 = (float *)malloc(sizeof(float) * tLen);

    roi1 = [img1 copy];
    [roi1 subImage:ref1];
    [roi1 crop:[roi1 xLoop] to:28 startAt:203];
    [roi1 crop:[roi1 yLoop] to:28 startAt:149];
    [roi1 saveAsKOImage:@"IMG_roi1.img"];

    roi2 = [img2 copy];
    [roi2 subImage:ref1];
    [roi2 crop:[roi2 xLoop] to:28 startAt:208];
    [roi2 crop:[roi2 yLoop] to:28 startAt:181];
    [roi2 saveAsKOImage:@"IMG_roi2.img"];

    n = [roi1 xDim] * [roi1 yDim];
    for (i = 0; i < tLen; i++) {
        p = [roi1 data] + i * n;
        sum1 = 0;
        for (j = 0; j < n; j++) {
            sum1 += p[j];
//    printf("%f\n", p[j]); // ???
        }
        ac1[i] = sum1;
//        printf("%d %f\n", i, sum1);
    }
    p = [roi2 data];
    for (i = 0; i < tLen; i++) {
        p = [roi2 data] + i * n;
        sum2 = 0;
        for (j = 0; j < n; j++) {
            sum2 += p[j];
        }
        ac2[i] = sum2;
    }

// accum
    for (i = 0; i < tLen; i++) {
        printf("%d %f %f\n", i, ac1[i], ac2[i]);
    }
    
    sum1 = sum2 = 0;
    for (i = 0; i < tLen; i++) {
        sum1 += ac1[i];
        sum2 += ac2[i];
        printf("%d %f %f\n", i, sum1/(i+1), sum2/(i+1));
    }
    free(ac1);
    free(ac2);
}

void
t1sim()
{
    int     iter, nIter = 4;
    int     i, n = 100;
    float   t, dt;
    float   TR = 3000;
//    float   TI = 2200;
//    float   Tread = 800;
    float   tcross = 0;
    float   T1 = 4000;  // water t1
    float   T2 = 3000;
    float   s;  // longitudinal mag

    s = 1.0;    // thermal eq
    dt = TR / n;
    for (iter = 0; iter < nIter; iter++) {
        tcross = 0;
        // inversion
        s = -s;
        for (i = 0; i < n; i++) {
            t = (float)i * TR / n;
            s = 1.0 - (1.0 - s) * exp(-dt / T1);
            if (tcross == 0 && s > 0) {
                tcross = t;
            }
            if (iter == (nIter - 1)) {
                printf("%f %f\n", t, s);
            }
        }
    }
    printf("tcross = %f\n", tcross);
}

void
partial()
{
    float   m1 = 0.5;
    float   m2 = 0.5;
    float   t2_1 = 60;
    float   t2_2 = 2000;
    float   te = 100;
    float   sg1, sg2;

    printf("partial\n");
    
    sg1 = m1 * exp(-te / t2_1); printf("1) %f %f\n", m1, sg1);
    sg2 = m2 * exp(-te / t2_2); printf("2) %f %f\n", m2, sg2);
    t2_1 = -te / log(sg1 + sg2); printf("t2_partial = %f\n", t2_1);
}
