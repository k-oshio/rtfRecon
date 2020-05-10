//
//
//

#import <NumKit/NumKit.h>
#import "RecImagePW.h"


typedef	struct gaEnt {
	float	th;
	int		origIx;
} gaEnt;

@implementation	RecImage (PW)

// Schofield paper
- (RecImage *)unwrapEst3D	// initial est using Laplacian
{
	RecImage	*est;
	RecImage	*cs, *sn, *ph1, *ph2;

	sn = [self copy];
	[sn thToSin];
	cs = [self copy];
	[cs thToCos];

	ph1 = [sn copy];
	[ph1 laplace3d:REC_FORWARD];
[ph1 saveAsKOImage:@"IMG_lp0.img"];
	[ph1 mulImage:cs];

	ph2 = [cs copy];
	[ph2 laplace3d:REC_FORWARD];
[ph2 saveAsKOImage:@"IMG_lp1.img"];
	[ph2 mulImage:sn];

	est = [ph1 copy];
	[est subImage:ph2];
	[est laplace3d:REC_INVERSE];

	return est;
}

- (void)phaseUnwrap3D
{
	RecImage	*mg, *phs, *mask;
	RecImage	*ph1, *ph3, *pnp;
	NSString	*path;
	int			i;

	mg = [self copy];
	[mg magnitude];
	phs = [self copy];
	[phs phase];

	mask = [self copy];
	[mask gauss3DLP:0.2];
	[mask magnitude];
	[mask thresEachSliceAt:0.1]; // slice by slice -> doesn't work very well
	[mask saveAsKOImage:@"IMG_mask.img"];

	[phs mulImage:mask];

	ph3 = [phs unwrapEst3D];

	[ph3 saveAsKOImage:@"IMG_est.img"];

exit(0);

//	[ph3 makeZeroMeanWithMask:mask];
//	[ph3 fermiWithRx:0.5 ry:0.5 d:0.1 x:0 y:0 invert:NO half:NO];	//##

	ph1 = [phs copy];	// initial
// iterative correction
// Schofield paper
	for (i = 0; i < 10; i++) {
		pnp = [ph3 copy];
		[pnp subImage:ph1];
//		[pnp makeZeroMeanWithMask:mask];
		path = [NSString stringWithFormat:@"IMG_np%d.img", i];
		[pnp saveAsKOImage:path];
		[pnp mod2PI];
//		path = [NSString stringWithFormat:@"IMG_mod%d.img", i];
		[pnp saveAsKOImage:path];
		[ph1 addImage:pnp];
		[ph1 makeZeroMeanWithMask:mask];
		path = [NSString stringWithFormat:@"IMG_ps%d.img", i];
		[ph1 saveAsKOImage:path];

//		ph3 = [ph1 unwrapEst];
	}
}

// make common interface with below
//- (RecLoop *)slcSft
- (RecLoop *)sftCenter:(RecLoop *)lp
{
	RecLoop		*newLp;
	float		pk;

	pk = [self findEchoCenterForLoop:lp];
	printf("sft:%f\n", pk);

	[self shift1d:lp by:(int)pk];     // npix shift

//	newLp = [self zeroFillToPo2:lp];
//    [self fft1d:newLp direction:REC_FORWARD];
    [self fft1d:lp direction:REC_FORWARD];

//	return newLp;
	return lp;
}

- (RecLoop *)sftCenter:(RecLoop *)lp zeroFillTo:(int)zDim
{
	RecLoop		*newLp;
	float		pk;

	pk = [self findEchoCenterForLoop:lp];
	printf("sft:%f\n", pk);

	[self shift1d:lp by:(int)pk];     // npix shift

	if (zDim > 0) {
		newLp = [self zeroFill:lp to:zDim];
		[self fft1d:newLp direction:REC_FORWARD];
		return newLp;
	} else {
		return lp;
	}
}

// ##### input is PW config
- (RecLoop *)halfFT2:(RecLoop *)lp
{
	RecLoop         *newLp;
	int             dim, newDim;
	RecImage        *phs, *neg_phs, *raw;
	RecLoopControl	*lc;
	int				i;
    BOOL            lower;
	int				st, len, hovr, ofs, hdim;
	float			pk;

	pk = [self findEchoCenterForLoop:lp];
	printf("sft:%f\n", pk);

	if (fabs(pk) < 4) {	// zero-fill and fft
		newLp = [self zeroFillToPo2:lp];
		[self shift1d:newLp by:(int)pk];     // npix shift
		[self fft1d:newLp direction:REC_FORWARD];
		return newLp;
	}
	// else half fourier
	dim = [lp dataLength];
	ofs = (int)pk;
	if (ofs > 0) {
		lower = YES;
		hovr = dim/2 - ofs;
	} else {
		lower = NO;
		hovr = dim/2 + ofs;
	}
	hdim = dim - hovr;
	newDim = hdim * 2;
	newDim = Rec_po2(newDim);
	newLp = [RecLoop loopWithDataLength:newDim];

	ofs = - (newDim/2 - hdim);
	[self replaceLoop:lp withLoop:newLp offset:ofs];
	[self setUnit:REC_FREQ forLoop:newLp];

    phs = [RecImage imageWithImage:self];
	lc = [self control];
	st = newDim/2 - hovr;
	len = hovr*2;
	[lc setRange:NSMakeRange(st, len) forLoop:newLp];
	[phs copyImage:self withControl:lc];
	// ft and take phase
	[phs fft1d:newLp direction:REC_FORWARD];	//###
	[phs toUnitImage];
	neg_phs = [phs copy];
	[neg_phs conjugate];

//=== iter loop ===
	raw = [self copy];	// initial
    if (lower) {
		st = 0;
		len = newDim/2 + hovr;
    } else {
		st = newDim/2 - hovr;
		len = newDim/2 + hovr;
    }
	for (i = 0; i < 2; i++) {
		// ft data
		[self fft1d:newLp direction:REC_FORWARD]; // forward
		[self multByImage:neg_phs];
		[self conjugate];
		[self multByImage:phs];	// re-apply phase
		[self fft1d:newLp direction:REC_INVERSE]; //
		// merge
		[lc resetRange];
		[lc setRange:NSMakeRange(st, len) forLoop:newLp];
		[self copyImage:raw withControl:lc];
	}
    [self fft1d:newLp direction:REC_FORWARD];

	return newLp;
}

// coil sensitivity map (complex)
- (RecImage *)ftFitWithMask:(RecImage *)mask
{
	RecImage	*mg, *est, *err;
	float		w = 0.02;
	int			i;

	// initial est
	mg = [self copy];
//	[mg magnitude];
	est = [mg copy];
	[est gauss3DLP:w];
	for (i = 0; i < 2; i++) {
	printf("iter = %d\n", i);
		err = [self copy];
		[err subImage:est];
		[err maskWithImage:mask];
		[err multByConst:0.2];
		[err gauss3DLP:w];
		[est addImage:err];
	}

	return est;
}

// img -> pw
- (RecImage *)reprojectWithMap:(RecImage *)radMap
{
    RecImage        *img, *pw, *proj;
    int             dim = [self dim];
    RecLoop         *pe, *rd;
    NSMutableArray  *loopArray;

//img  [ch, sl, ky, kx]
//pw   [ch, pe, sl, rd]
//proj [ch, sl, pe, rd]
//radMap:      [pe, rd]

	img = [self copy];
    [img makeComplex];	
    [img fft2d:REC_INVERSE];      // img -> raw

    loopArray = [NSMutableArray arrayWithArray:[self loops]];   // [ch, sl, ky, kx]
    pe = [radMap yLoop];
    rd = [radMap xLoop];
    [loopArray replaceObjectAtIndex:(dim - 2) withObject:pe];
    [loopArray replaceObjectAtIndex:(dim - 1) withObject:rd];
                                // [ch, sl, pe, rd]
    proj = [RecImage imageOfType:RECIMAGE_COMPLEX withLoopArray:loopArray];
    [proj copyUnitOfImage:img];
    [proj setUnit:REC_FREQ forLoop:pe];
    [proj setUnit:REC_FREQ forLoop:rd];
    [proj resample:img withMap:radMap];
    [proj fft1d:rd direction:REC_FORWARD];

    [loopArray exchangeObjectAtIndex:(dim - 3) withObjectAtIndex:(dim - 2)];
                                // [ch, pe, sl, rd]
    pw = [RecImage imageOfType:RECIMAGE_COMPLEX withLoopArray:loopArray];
    [pw copyImage:proj];			// single-ch, complex

    return pw;
}

- (RecImage *)pDiffWith:(RecImage *)ref
{
	RecImage	*img1, *img2;

// avoid side effect
	img1 = [self copy];
	img2 = [ref copy];

	[img2 conjugate];
	[img1 multByImage:img2];
	
	return img1;
}

#define NT 8
// === plans ===
// *2-D phase error for ref views are now OK
// use 0th and 1st phase + trig interp for first implementation
// implement self-ref + 2-D correction later
//		chk corresponding view # for g.a. -> view table ?
//
- (void)radPhaseCorrWithGATab:(RecImage *)tab  // in/out : raw [ch sl proj rd]
{
	RecImage		*img[NT];
	RecImage		*corr;
//	NSPoint			pt;

	RecLoop			*xLoop, *yLoop, *zLoop, *ch;
	RecLoopControl	*lc;		// full
	int				i, len, n, view, nView;
	float			*re, *im;
	float			*phs0_p, *phs1_p;
	float			*gaTab = [tab data];
	RecImage		*phs0, *phs1, *phs0Low, *phs1Low;
	RecLoopIndex	*viewIndex;

	[self removePointLoops];	// Echo
	xLoop = [self xLoop];
	yLoop = [self yLoop];
	zLoop = [self zLoop];
	nView = [yLoop dataLength];

	for (i = 0; i < NT; i++) {
		img[i] = [self sliceAtIndex:nView - NT + i forLoop:yLoop];
		[img[i] fGauss2DLP:0.1];
		if (i >= NT/2) {
			[img[i] xFlip];
		}
	}

//	ch = [RecLoop findLoop:@"Channel"];	// GE version
	ch = [self topLoop];	// loop ID should be common among manufacturers
	if (ch == [self zLoop]) {
		printf("single channel data\n");
		exit(0);	// implement single ch later
	}

	for (i = 0; i < 8; i++) {
		[img[i] fft2d:REC_FORWARD];
	}

	xLoop = [self xLoop];
	len = [self xDim];
	phs0Low = [RecImage imageOfType:RECIMAGE_REAL xDim:NT+1];
	phs1Low = [RecImage imageOfType:RECIMAGE_REAL xDim:NT+1];
	phs0_p = [phs0Low data];
	phs1_p = [phs1Low data];

	for (i = 0; i < NT/2; i++) {
		NSString *path;
		float		*p, *q;
		float		a0, a1;

		corr = [img[i+NT/2] pDiffWith:img[i]];
		corr = [corr avgForLoop:ch];
		corr = [corr avgForLoop:zLoop];
	//	[corr fft1d:xLoop direction:REC_INVERSE];
		path = [NSString stringWithFormat:@"corr%d.img", i];
	//	[corr  saveAsKOImage:path];

		p = [corr data];
		q = p + [corr dataLength];
		a1 =  Rec_est_1st(p, q, 1, len);
		Rec_corr1st(p, q, 1, len, a1/2);		// this is necessary for 0th est
		a0 =  Rec_est_0th(p, q, 1, len);
//		Rec_corr0th(p, q, 1, len, a0);		// this is not necessary for the result
//	printf("%d %f %f\n", i, a0, a1);
//		path = [NSString stringWithFormat:@"corr_corr%d.img", i];
//		[corr  saveAsKOImage:path];

		phs0_p[i] = -a0/2;
		phs1_p[i] = -a1/2;
		phs0_p[i+NT/2] = -phs0_p[i];
		phs1_p[i+NT/2] = -phs1_p[i];

	}
//	for (i = 0; i < NT; i++) {
//		printf("%d %f %f\n", i * nView/NT, phs0_p[i], phs1_p[i]);
//	}
	phs0_p[NT] = phs0_p[0];

// === correction
	phs0 = [phs0Low scale1dLoop:[phs0Low xLoop] by:(float)nView/NT to:len];
	phs0_p = [phs0 data];
	phs1 = [phs1Low scale1dLoop:[phs1Low xLoop] by:(float)nView/NT to:len];
	phs1_p = [phs1 data];

//[phs0 saveAsKOImage:@"phs0.img"];
//[phs1 saveAsKOImage:@"phs1.img"];

	[self fft1d:xLoop direction:REC_FORWARD];
    lc = [self outerLoopControl];
	viewIndex = [lc loopIndexForLoop:yLoop];
    [lc rewind];
    n = [lc loopLength];
    for (i = 0; i < n; i++) {
		view = [viewIndex current];	// look up yLoop
		view = (int)(gaTab[view] * nView / M_PI / 2);
//	printf("%d\n", view);
        re = [self currentDataWithControl:lc];
        im = re + [self dataLength];
        Rec_corr0th(re, im, 1, len, phs0_p[view]);	// center: (N-1)/2, phs: radian / FOV
        Rec_corr1st(re, im, 1, len, phs1_p[view]);	// center: (N-1)/2, phs: radian / FOV

        [lc increment];
    }
	[self fft1d:xLoop direction:REC_INVERSE];	// sin -> raw
}

// NT -> yDim
// 0th order not correct #### 12-15-2017
- (void)radPhaseCorr360	// self-pcorr using 360 data (g.a. or not), input is raw
{
	RecImage		*img1, *img2, *corr;
//	NSString		*path;
	RecLoop			*xLoop, *yLoop, *zLoop, *ch;
	RecLoopControl	*lc;		// full
	int				i, len, n, view, nView;
	float			*re, *im;
	float			*phs0_p, *phs1_p;
	RecImage		*phs0, *phs1;
	RecLoopIndex	*viewIndex;

	[self removePointLoops];	// Echo
	xLoop = [self xLoop];
	yLoop = [self yLoop];
	zLoop = [self zLoop];
	nView = [yLoop dataLength];

	ch = [self topLoop];	// loop ID should be common among manufacturers
	if (ch == [self zLoop]) {
		printf("single channel data\n");
		exit(0);	// implement single ch later
	}

// 0-th correction & 1-st estimation
// (error due to x-shift... should not exist)
	xLoop = [self xLoop];
	len = [self xDim];
	nView = [self yDim];
	phs0 = [RecImage imageOfType:RECIMAGE_REAL xDim:nView];
	phs1 = [RecImage imageOfType:RECIMAGE_REAL xDim:nView];
	phs0_p = [phs0 data];
	phs1_p = [phs1 data];

	for (i = 0; i < nView/2; i++) {
		float		*p, *q;
		float		a0, a1;

		img1 = [self sliceAtIndex:i forLoop:yLoop];
		img2 = [self sliceAtIndex:i + nView/2 forLoop:yLoop];
		// flip img2 (k-space)
		[img2 xFlip];

		// pDiff (sinogram)
		[img1 fft1d:xLoop direction:REC_FORWARD];
		[img2 fft1d:xLoop direction:REC_FORWARD];

		corr = [img1 pDiffWith:img2];
		corr = [corr avgForLoop:ch];
		corr = [corr avgForLoop:zLoop];

		p = [corr real];
		q = [corr imag];
		a1 =  Rec_est_1st(p, q, 1, len) / 2;
		Rec_corr1st(p, q, 1, len, a1); //a1/2);		// this is necessary for 0th est
		a0 =  Rec_est_0th(p, q, 1, len);
	//	Rec_corr0th(p, q, 1, len, a0/2);		// this is not necessary for the result

		phs0_p[i] = a0;
		phs1_p[i] = a1; //a1;
		phs0_p[i+nView/2] = -a0;
		phs1_p[i+nView/2] = a1; //a1;
//printf("%d %f %f\n", i, a1, a0);
	}
	// divide phase by 2 (1D unwrap, then div)
	Rec_unwrap_1d(phs0_p, nView, 1);
	for (i = 0; i < nView; i++) {
		phs0_p[i] /= 2.0;
	}

	// correction
	[self fft1d:xLoop direction:REC_FORWARD];
    lc = [self outerLoopControl];
	viewIndex = [lc loopIndexForLoop:yLoop];
    [lc rewind];
    n = [lc loopLength];
    for (i = 0; i < n; i++) {
		view = [viewIndex current];	// look up yLoop
        re = [self currentDataWithControl:lc];
        im = re + [self dataLength];
        Rec_corr0th(re, im, 1, len, phs0_p[view]);	// center: (N-1)/2, phs: radian / FOV
        Rec_corr1st(re, im, 1, len, phs1_p[view]);	// center: (N-1)/2, phs: radian / FOV

        [lc increment];
    }
	[self fft1d:xLoop direction:REC_INVERSE];	// sin -> raw
}

// clockwise ...
float toshibaRad(float th)
{
	return -0.5 * M_PI - th;
}

- (void)initRadialTrajWithTab:(RecImage *)tab startPoint:(int)st
{
	int				nproj, xdim;
	RecLoopControl	*lc;
	int				i, j;
	float			*tabp;
	float			r, th, cs, sn;
	float			kx, ky, den;
	float			*pkx, *pky, *pden, *pwt;

	xdim  = [self xDim];
	nproj = [self yDim];
	tabp = [tab data];

	lc = [self outerLoopControl];
	[lc rewind];

    for (i = 0; i < nproj; i++) {
		th = tabp[i];
	//	th = toshibaRad(th);	// ##
        cs = cos(th);
        sn = sin(th);
		pkx = [self currentDataWithControl:lc];
		pky = pkx + dataLength;
		pden = pky + dataLength;
		pwt = pden + dataLength;
        for (j = 0; j < xdim; j++) {
            r = (j - xdim/2.0) / xdim;	// center is N/2
            kx = r * cs;
            ky = r * sn;
            den = fabs(r) * xdim / nproj;
            if (den < 1.0 / nproj) {
                den = 1.0 / nproj;
            }
			if (j < st) {
				den = 0;
			}
            pkx[j] = kx;        // [-0.5..0.5]
            pky[j] = ky;        // [-0.5..0.5]
            pden[j] = den;		// initial value... replaced by iterative density est
			pwt[j] = 1.0;
        }
		[lc increment];
    }
}

- (void)reorder:(RecLoop *)lp withTab:(RecImage *)tab
{
	int			n = [lp dataLength];
	float		*buf;
	float		*origIx = [tab imag];
    void		(^proc)(float *p, int len, int skip);

	buf = (float *)malloc(sizeof(float) * n);

    proc = ^void(float *p, int len, int skip) {
        int     i, ix;
		
		for (i = ix = 0; i < len; i++, ix += skip) {
			buf[i] = p[ix];
		}
		for (i = ix = 0; i < len; i++, ix += skip) {
			p[ix] = buf[(int)(origIx[i])];
		//	p[ix] = buf[i];
		}
    };
    [self apply1dProc:proc forLoop:lp];
	free(buf);
}

- (void)invReorder:(RecLoop *)lp withTab:(RecImage *)tab
{
	int			n = [lp dataLength];
	float		*buf;
	float		*origIx = [tab imag];
    void		(^proc)(float *p, int len, int skip);

	buf = (float *)malloc(sizeof(float) * n);

    proc = ^void(float *p, int len, int skip) {
        int     i, ix;
		
		for (i = ix = 0; i < len; i++, ix += skip) {
			buf[(int)(origIx[i])] = p[ix];
		}
		for (i = ix = 0; i < len; i++, ix += skip) {
			p[ix] = buf[i];
		}
    };
    [self apply1dProc:proc forLoop:lp];
	free(buf);
}

- (void)initRadialTrajWithSortTab:(RecImage *)tab startPoint:(int)st
{
	int				nproj, xdim;
	RecLoopControl	*lc;
	int				i, j;
	float			*tabp;
	float			r, th, cs, sn;
	float			kx, ky, den;
	float			*pkx, *pky, *pden, *pwt;

	xdim  = [self xDim];
	nproj = [self yDim];
	tabp = [tab data];

	lc = [self outerLoopControl];
	[lc rewind];

    for (i = 0; i < nproj; i++) {
		th = tabp[i];	// only diff with initRadialTrajWithTab
	//	th = toshibaRad(th);	// ##
        cs = cos(th);
        sn = sin(th);
		pkx = [self currentDataWithControl:lc];
		pky = pkx + dataLength;
		pden = pky + dataLength;
		pwt = pden + dataLength;
        for (j = 0; j < xdim; j++) {
            r = (j - xdim/2.0) / xdim;	// center is N/2
            kx = r * cs;
            ky = r * sn;
            den = fabs(r) * xdim / nproj;
            if (den < 1.0 / nproj) {
                den = 1.0 / nproj;
            }
			if (j < st) {
				den = 0;
			}
            pkx[j] = kx;        // [-0.5..0.5]
            pky[j] = ky;        // [-0.5..0.5]
            pden[j] = den;		// initial value... replaced by iterative density est
			pwt[j] = 1.0;
        }
		[lc increment];
    }
}

- (void)initRadialMapWithTab:(RecImage *)tab	// always full echo
{
	int				nproj, xdim;
	RecLoopControl	*lc;
	int				i, j;
	float			*tabp;
	float			r, th, cs, sn;
	float			kx, ky;
	float			*pkx, *pky;

	xdim  = [self xDim];
	nproj = [self yDim];
	tabp = [tab data];

	lc = [self outerLoopControl];
	[lc rewind];

    for (i = 0; i < nproj; i++) {
		th = tabp[i];
		th = toshibaRad(th);	// ##
        cs = cos(th);
        sn = sin(th);
		pkx = [self currentDataWithControl:lc];
		pky = pkx + dataLength;
        for (j = 0; j < xdim; j++) {
            r = (j - xdim/2.0) / xdim;	// center is N/2
            kx = r * cs;
            ky = r * sn;
            pkx[j] = kx;        // [-0.5..0.5]
            pky[j] = ky;        // [-0.5..0.5]
        }
		[lc increment];
    }
}

// pw -> sinogram
- (RecImage *)pwToSin
{
    RecImage        *raw;
    int             dim = [self dim];
    NSMutableArray  *loopArray;

    loopArray = [NSMutableArray arrayWithArray:[self loops]];   // [ch, pe, sl, rd]
    [loopArray exchangeObjectAtIndex:(dim - 3) withObjectAtIndex:(dim - 2)];
    raw = [RecImage imageOfType:RECIMAGE_COMPLEX withLoopArray:loopArray];
    [raw copyImage:self];

    return raw;
}

// ### document traj and map data -> merge or auto-conversion
// radial (y = theta)
- (RecImage *)mapForRadial
{
	RecImage	*map;
    int         xDim = [self xDim];
    int         yDim = [self yDim];
	float		*mapX, *mapY;
	int			i, j, ix, len;
	float		r, th;
    float       cs, sn;

//	map = [RecImage imageOfType:RECIMAGE_MAP xDim:xDim yDim:yDim];	// outer loop first
	map = [RecImage imageWithImage:self];
	len = [map dataLength];
	mapX = [map data];
	mapY = mapX + len;
	
	for (i = ix = 0; i < yDim; i++) {
		th = i * M_PI / yDim;
        cs = cos(th);
        sn = sin(th);
		for (j = 0; j < xDim; j++, ix++) {
			r = ((float)j - xDim/2) / xDim;	//	[-0.5 .. 0.5]
			mapX[ix] = r * cs;
			mapY[ix] = r * sn;
		}
	}
	return map;
}



- (void)toPolarFrom:(RecImage *)img
{
    RecImage        *map = [self mapForPolarWithRMin:1 logR:NO];
//    RecImage        *map = [self mapForPolar];
    [self resample:img withMap:map];
}

- (RecImage *)toCoronal
{
	RecImage	*cor;
	cor = [RecImage imageWithImage:self];
	[cor swapLoop:[self zLoop] withLoop:[self yLoop]];
	[cor copyImage:self];
	return cor;
}

- (RecImage *)toSagittal
{
	RecImage	*sag;
	sag = [RecImage imageWithImage:self];
	[sag swapLoop:[self zLoop] withLoop:[self yLoop]];
	[sag swapLoop:[self xLoop] withLoop:[self yLoop]];
	[sag copyImage:self];
	return sag;
}

- (void)pwWinWithRx:(float)rx ry:(float)ry d:(float)d x:(float)xc y:(float)yc
{
	RecImage	*win;
	RecLoop		*zLp = [self zLoop];
	RecLoop		*xLp = [self xLoop];
	int			zDim = [self zDim];
	int			xDim = [self xDim];
	float		*p;
	int			i, j;
	float		x, th, x0, r, w;
	float		cs, sn, cp, cm;

	win = [RecImage imageOfType:RECIMAGE_REAL withLoops:zLp, xLp, nil];
	p = [win data];

	for (i = 0; i < zDim; i++) {
		th = (float)i * M_PI / zDim;
		cs = cos(th); sn = sin(th);
		cp = cs * cs;
		cm = 1 - cp;
		x0 = xc * cs - yc * sn;
		r = rx * cp + ry * cm; 
		for (j = 0; j < xDim; j++) {
			x = ((float)j - xDim/2) * 2 / xDim - x0;
			x = fabs(x) - r;
			if (x < 0) {
				w = 1.0;
			} else
			if (x > d) {
				w = 0;
			} else {
				w = cos(x * M_PI / d) * 0.5 + 0.5;
			}
			p[i * xDim + j] = w;
		}
	}
	[win saveAsKOImage:@"pwwin.img"];
	[self multByImage:win];
}

- (void)imgWin
{
	[self fGauss2DLP:0.25];
//	[self fermiWithW:0.35 andD:0.3];
}

- (void)SATx:(float)xofs w:(float)w th:(float)th
{
	float		rr, ww;
	float		x, y;
	int			i, j, xDim, yDim;
	RecLoop		*xLp, *yLp;
	RecImage	*win;
	float		*p;

	xLp = [self xLoop];
	yLp = [self yLoop];
	win = [RecImage imageOfType:RECIMAGE_REAL withLoops:yLp, xLp, nil];
	p = [win data];
	xDim = [win xDim];
	yDim = [win yDim];

	// set win
	th *= -M_PI / 180.0;
	for (i = 0; i < yDim; i++) {
		y = ((float)i - yDim/2);
		for (j = 0; j < xDim; j++) {
			x = ((float)j - xDim/2);
			rr = x * cos(th) - y * sin(th) - xofs;
			if (fabs(rr) < w/2) {
				ww = 0;
			} else {
				ww = 1.0;
			}
			p[i * xDim + j] = ww;
		}
	}
	[self multByImage:win];
}

// unit of shift: pixel
// self = corr image [ch, proj, z, x]
// output: [proj] (z-shift only)
- (RecImage *)corrToSft	// top level
{
	RecImage	*sft;

	sft = [self corrToSftMC];
	return [sft combineSftXY];
}

- (RecImage *)corrToSftZ	// top level
{
	RecImage	*sft;

	sft = [self corrToSftMC];
	return [sft combineSftZ];
}

- (RecImage *)corrToSftMC	// multi-channel, xy
{
    RecImage        *sft, *cVal;
	RecLoopControl	*lc;
	RecLoop			*tLoop;
	int				xDim, yDim, tDim, i, nImages;
	NSPoint			pt;
	float			*p, *x, *y, *m, mx;

	tLoop = [self zLoop];
	lc = [self control];
	[lc deactivateXY];
	sft  = [RecImage imageOfType:RECIMAGE_MAP withControl:lc];	// peak position
	cVal = [RecImage imageOfType:RECIMAGE_REAL withControl:lc];	// peak correlation value

	x = [sft real];
	y = [sft imag];
	m = [cVal real];
	xDim = [self xDim];
	yDim = [self yDim];
	tDim = [self zDim];
	nImages = [self nImages];	// tDim * nCh
//	lc = [self control];
//	[lc deactivateXY];

	for (i = 0; i < nImages; i++) {
		p = [self currentDataWithControl:lc];
		// unit: pixels
		pt = Rec_find_peak2_mx(p, xDim, yDim, &mx);
		// remove erroneous val
		// this step is necessary
		if (fabs(pt.x) > 4.0) {
			pt.x = 0;
		}
		if (fabs(pt.y) > 4.0) {
			pt.y = 0;
		}
		x[i] = pt.x;
		y[i] = pt.y;
		m[i] = mx;
//	printf("%d %f %f %f\n", i, pt.x, pt.y, mx);
		[lc increment];
	}

    return sft;
}

- (RecImage *)combineSftXY	// svd or weighed sum		### not implemented yet
{
	RecImage		*sft;
	RecImage		*U, *S, *Vt;			// SVD
	RecLoop			*tLoop;
	int				tDim, chDim, i;
	float			*p, *q, *x, *y, mx;

	sft = [self copy];
//	[sft takeImagPart];		// z-shift only ### complex svd not working
	[sft trans];

	[sft svd_U:&U S:&S Vt:&Vt];	// dbg svd (mem alloc issue) ###
	[U trans];
	[U saveAsKOImage:@"IMG_U.img"];
	[S saveAsKOImage:@"IMG_s.img"];
	[Vt saveAsKOImage:@"IMG_Vt.img"];

//	copy U to sft
	tLoop = [self xLoop];
	chDim = [self yDim];
	tDim = [tLoop dataLength];
	sft = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:tLoop, nil];
	x = [sft real];
	y = [sft imag];
	p = [U real];
	q = [U imag];
	mx = 0;
	for (i = 0; i < chDim; i++) {
		mx += [Vt data][i];
	}
	mx /= chDim;
	mx *= [S data][0];
	for (i = 0; i < tDim; i++) {
		x[i] = mx * p[i];
		y[i] = mx * q[i];
	//	printf("%d %f\n", i, y[i]);
	}
[sft saveAsKOImage:@"sft_xy.img"];

    return sft;
}

// self is sft
- (RecImage *)combineSftZ;	// svd ? imag only			### not implemented yet
{
	RecImage		*sft;
	RecImage		*U, *S, *Vt;			// SVD
	RecLoop			*tLoop;
	int				tDim, chDim, i;
	float			*p, *y, mx;

	sft = [self copy];
	[sft takeImagPart];		// z-shift only ### complex svd not working
	[sft trans];

	[sft svd_U:&U S:&S Vt:&Vt];	// dbg svd (mem alloc issue) ###
	[U trans];
	[U saveAsKOImage:@"IMG_U.img"];
	[S saveAsKOImage:@"IMG_s.img"];
	[Vt saveAsKOImage:@"IMG_Vt.img"];

//	copy U to sft
	tLoop = [self xLoop];
	chDim = [self yDim];
	tDim = [tLoop dataLength];
	sft = [RecImage imageOfType:RECIMAGE_REAL withLoops:tLoop, nil];
	y = [sft data];
	p = [U data];
	mx = 0;
	for (i = 0; i < chDim; i++) {
		mx += [Vt data][i];
	}
	mx /= chDim;
	mx *= [S data][0];
	for (i = 0; i < tDim; i++) {
		y[i] = mx * p[i];
	//	printf("%d %f\n", i, y[i]);
	}
//[sft saveAsKOImage:@"sft_svd.img"];

    return sft;
}

- (RecImage *)detectShift			// optical flow or laplacian
{
	RecImage	*sft, *nav;
	int			i, j, xDim, yDim;
	float		*p, *sf, mn;

	nav = [self copy];
	[nav magnitude];	// self is complex
// baseline removal along y (implement cos filter later)
//	[nav cosFilter:[nav yLoop] order:2 keepDC:YES];
	[nav gauss1DLP:0.1 forLoop:[nav xLoop]];
	[nav normalizeForLoop:[nav xLoop]];

	xDim = [nav xDim];
	yDim = [nav yDim];
//	ref = [nav avgForLoop:[nav yLoop]];
[nav saveAsKOImage:@"nav.img"];
//[ref saveAsKOImage:@"ref.img"];

	[nav laplace1dForLoop:[nav xLoop] direction:REC_FORWARD];

	// select lower half
	p = [nav data];
	for (i = 0; i < yDim; i++) {
		for (j = 0; j < xDim; j++) {
			if (j < xDim/2) {
				p[i * xDim + j] = 0;
			}
		}
	}
[nav saveAsKOImage:@"nav_lap.img"];

	sft = [RecImage imageOfType:RECIMAGE_REAL withLoops:[nav yLoop], nil];
	sf = [sft data];

	for (i = 0; i < yDim; i++) {
		sf[i] = Rec_find_peak(p + i * xDim, 1, xDim);
	//	printf("%d %f\n", i, sf[i]);
	}
	mn = [sft meanVal];
	[sft addConst:-mn];

[sft saveAsKOImage:@"sft.img"];

	return sft;
}

// 1-D optical flow
- (RecImage *)opticalFlow1d
{
	RecImage	*sft, *nav, *ref, *var, *dif;
	int			i, j, xDim, yDim;
	int			ix, n;
	float		*p, *q, *gr, *sf, thres;

	nav = [self copy];
	[nav magnitude];	// self is complex
// baseline removal along y (implement cos filter later)
//	[nav cosFilter:[nav yLoop] order:2 keepDC:YES];
	[nav gauss1DLP:0.1 forLoop:[nav xLoop]];
	[nav normalizeForLoop:[nav xLoop]];

	xDim = [nav xDim];
	yDim = [nav yDim];
	ref = [nav avgForLoop:[nav yLoop]];
[nav saveAsKOImage:@"nav.img"];
[ref saveAsKOImage:@"ref.img"];

// determine window center and make gaussian window
//	vr = [nav varForLoop:[nav yLoop]];
	dif = [RecImage imageWithImage:nav];
	p = [nav data];
	q = [dif data];
	for (i = 0; i < yDim; i++) {
		for (j = 1; j < xDim; j++) {
			ix = i * xDim + j;
			q[ix] = p[ix] - p[ix - 1];
		}
	}
//	[dif normalizeForLoop:[dif xLoop]];
[dif saveAsKOImage:@"df.img"];

// calc optical flow with ref (avg)
	sft = [RecImage imageOfType:RECIMAGE_REAL withLoops:[nav yLoop], nil];
	sf = [sft data];
	p = [nav data];
	q = [ref data];
	gr = [dif data];
	thres = [dif maxVal] * 0.6;
	for (i = 0; i < yDim; i++) {
		n = 0;
		sf[i] = 0;
		for (j = 0; j < xDim; j++) {
			ix = i * xDim + j;
			if (gr[ix] < thres) continue;
			sf[i] += (p[ix] - q[j]) / gr[ix];
			n++;
		}
		sf[i] /= n;
	//	printf("%d %f %d\n", i, sf[i], n);
	}
	[sft cosFilter:[sft xLoop] order:2 keepDC:NO];
[sft saveAsKOImage:@"sft.img"];

	return sft;
}

// not done yet ##
- (RecImage *)opticalFlow2dWithRef:(RecImage *)ref
{
	RecImage	*img;

	return img;
}

// not done yet ##
- (RecImage *)opticalFlow3dWithRef:(RecImage *)ref
{
	RecImage	*v;
	RecImage	*gx, *gy, *gz, *gt;
	RecImage	*xx, *yy, *zz, *xy, *xz, *yz;
	RecImage	*bx, *by, *bz;
	Num_mat		*M, *b;

	int			i, n;
	int			sts, nerr;
	float		*pxx, *pyy, *pzz, *pxy, *pxz, *pyz;
	float		*pbx, *pby, *pbz;
	float		*pvx, *pvy, *pvz;
	float		w = 0.03;

// gradient
	[ref copyLoopsOf:self];
	gx = [ref copy];
	[gx grad1dForLoop:[gx xLoop]];
	gy = [ref copy];
	[gy grad1dForLoop:[gy yLoop]];
	gz = [ref copy];
	[gz grad1dForLoop:[gz zLoop]];
	gt = [self copy];
	[gt subImage:ref];

//	[gx saveAsKOImage:@"img_gx.img"];
//	[gy saveAsKOImage:@"img_gy.img"];
//	[gz saveAsKOImage:@"img_gz.img"];
//	[gt saveAsKOImage:@"img_dt.img"];

// M components
	xx = [gx copy];
	[xx mulImage:gx];
	yy = [gy copy];
	[yy mulImage:gy];
	zz = [gz copy];
	[zz mulImage:gz];
	xy = [gx copy];
	[xy mulImage:gy];
	xz = [gx copy];
	[xz mulImage:gz];
	yz = [gy copy];
	[yz mulImage:gz];

	[xx gauss3DLP:w];
	[yy gauss3DLP:w];
	[zz gauss3DLP:w];
	[xy gauss3DLP:w];
	[xz gauss3DLP:w];
	[yz gauss3DLP:w];

//	[xx saveAsKOImage:@"img_xx.img"];
//	[yy saveAsKOImage:@"img_yy.img"];
//	[zz saveAsKOImage:@"img_zz.img"];
//	[xy saveAsKOImage:@"img_xy.img"];
//	[xz saveAsKOImage:@"img_xz.img"];
//	[yz saveAsKOImage:@"img_yz.img"];


// b components
	bx = [gt copy];
	[bx mulImage:gx];
	by = [gt copy];
	[by mulImage:gy];
	bz = [gt copy];
	[bz mulImage:gz];

	[bx gauss3DLP:w];
	[by gauss3DLP:w];
	[bz gauss3DLP:w];

//	[bx saveAsKOImage:@"img_bx.img"];
//	[by saveAsKOImage:@"img_by.img"];
//	[bz saveAsKOImage:@"img_bz.img"];


// v image
	v = [RecImage imageOfType:RECIMAGE_VECTOR withImage:self];

// copy to matrix, invert, copy back to image
	M = Num_new_mat(3, 3);
	b = Num_new_mat(3, 1);

	pxx = [xx data];
	pyy = [yy data];
	pzz = [zz data];
	pxy = [xy data];
	pxz = [xz data];
	pyz = [yz data];

	pbx = [bx data];
	pby = [by data];
	pbz = [bz data];

	n = [self dataLength];
	pvx = [v data];
	pvy = pvx + n;
	pvz = pvy + n;

	nerr = 0;
	for (i = 0; i < n; i++) {
		float	mx = 20;
	// copy to mat
		M->data[0]				= pxx[i];
		M->data[4]				= pyy[i];
		M->data[8]				= pzz[i];
		M->data[1] = M->data[3] = pxy[i];
		M->data[2] = M->data[6] = pxz[i];
		M->data[5] = M->data[7] = pyz[i];
		b->data[0] = pbx[i];
		b->data[1] = pby[i];
		b->data[2] = pbz[i];

		// LAPACK
		sts = Num_inv(M, b, b);	// not tested ###

		if (sts != 0 ||
			fabs(b->data[0]) > mx ||
			fabs(b->data[1]) > mx ||
			fabs(b->data[2]) > mx ) {
			pvx[i] = pvy[i] = pvz[i] = 0;
			nerr++;
		} else {
	// copy to v
			pvx[i] = b->data[0];
			pvy[i] = b->data[1];
			pvz[i] = b->data[2];
		}
	}

	Num_free_mat(M);
	Num_free_mat(b);

//	[v saveAsKOImage:@"img_v.img"];
//	printf("sing: %f\n", (float)nerr/n);
	printf("v range: %f / %f\n", [v maxVal], [v minVal]);

	return v;
}

// shift estimation from center-of-kspace view (POCS)
// input is PW (ch, z:view, y:slc, x:read)
// output unit is pixels
- (RecImage *)shiftFromK0
{
	RecImage	*prj;
	RecImage	*mean;
	RecImage	*st, *mv, *sft, *nsft;		// stationary[xDim], moving[xDim], shift[yDim]
	RecImage	*est, *dif;
	int			i, j, ix, iter, nIter = 20; //20;
	int			xDim, yDim;
	float		*m, *p;
	float		err, mx;
	NSString	*path;
	BOOL		dbg = YES;

// ============ make 2D input from PW ======
	prj = [self combineForLoop:[self topLoop]];	// PW
	prj = [prj sumForLoop:[prj xLoop]];
	[prj removePointLoops];
	[prj takeRealPart];

// ============ main proc start ======
	[prj dumpLoops];
	xDim = [prj xDim];
	yDim = [prj yDim];

// normalize (mag)
	mean = [prj avgForLoop:[prj xLoop]];
	m = [mean data];
	p = [prj data];
	for (i = 0; i < yDim; i++) {
		for (j = 0; j < xDim; j++) {
			if (m[i] != 0) {
				p[i * xDim + j] /= m[i];
			}
		}
	}

// gradient along x
	dif = [RecImage imageWithImage:prj];
	p = [prj data];
	m = [dif data];
	for (i = 0; i < yDim; i++) {
		for (j = 1; j < xDim; j++) {
			ix = i * xDim + j;
			m[ix] = p[ix] - p[ix - 1];
		}
	}
	prj = [dif copy];

// tri win
	[prj fTriWin1DforLoop:[prj xLoop]];
    
// remove low freq along y
    [prj cosFilter:[prj yLoop] order:12 keepDC:YES]; // 4

// input
    if (dbg) {
        [prj saveAsKOImage:@"IMG_in"];
    }

// initial est
    int init_mode = 2;
    switch (init_mode) {
    case 0 :
        mv = [prj sliceAtIndex:0 forLoop:[prj yLoop]];
        st = [RecImage imageWithImage:mv];  // 0
        break;
    case 1 :
        st = [[prj avgForLoop:[prj yLoop]] multByConst:0.5];  // p/2
        mv = [st copy];
        break;
    case 2 :
        st = [[prj avgForLoop:[prj yLoop]] multByConst:0.9];  // p avg
        mv = [prj sliceAtIndex:0 forLoop:[prj yLoop]];
        [mv subImage:st];
        break;
    case 3 :
        st = [prj avgForLoop:[prj yLoop]];  // p avg
        mv = [[prj subImage:st] avgForLoop:[prj yLoop]];
        break;
    }

    // initial
    if (dbg) {
        [st saveAsKOImage:@"IMG_st0"];
        [mv saveAsKOImage:@"IMG_mv0"];
        [prj subImage:st];
        [prj saveAsKOImage:@"IMG_in_st"];
    }

// === iteration ===
	for (iter = 1; iter <= nIter; iter++) {
        int lp_mode = 0;

		est = [prj copy];
		[est subImage:st];
        
		sft = [est corr1dWithRef:mv];	// unit:pixels
		mx = [sft meanVal];
		[sft addConst:-mx];             // zero-mean

        nsft = [sft copy];
		[nsft negate];
		
		[est clear];
		[est copyImage:mv];
		est = [est correctShift1d:nsft forLoop:[mv xLoop]]; // pixels

		[est addImage:st];
		dif = [prj copy];
		[dif subImage:est];
		[dif fTriWin1DforLoop:[dif xLoop]];

		err = [dif rmsVal];
        if (dbg) {
            printf("%d %e\n", iter, err);
        }
		mean = [dif avgForLoop:[dif yLoop]];
	    [st addImage:mean]; // ok

		mv = [prj copy];
	 	[mv subImage:st];
		mv = [mv correctShift1d:sft forLoop:[mv xLoop]];
		mv = [mv avgForLoop:[mv yLoop]];

	//    if (dbg && (iter == nIter)) {
            if (dbg) {
                ix = 1;
            path = [NSString stringWithFormat:@"IMG_sft%d.img", ix];
            [sft saveAsKOImage:path];
			path = [NSString stringWithFormat:@"IMG_mv%d.img", ix];
			[mv saveAsKOImage:path];
			path = [NSString stringWithFormat:@"IMG_st%d.img", ix];
			[st saveAsKOImage:path];
			path = [NSString stringWithFormat:@"IMG_dif%d.img", ix];
			[dif saveAsKOImage:path];
			path = [NSString stringWithFormat:@"IMG_est%d.img", ix];
			[est saveAsKOImage:path];
            dif = [est copy];
            [dif subImage:st];
            path = [NSString stringWithFormat:@"IMG_ms%d.img", ix];
            [dif saveAsKOImage:path];

            p = [sft data];
            for (i = 0; i < [sft xDim]; i++) {
                printf("%d %f\n", i, p[i]);
            }
		}
	}

	return sft;
}

- (RecImage *)stepCorrWithPW:(RecImage *)pw gridder:(RecGridder *)grid sft:(RecImage *)sft // input(self) is img
{
    RecImage        *img, *imgs, *mxSft;
    RecImage        *sft_scl, *pws, *raw;
    RecImage        *tmp_img, *img_c;
    int             i;
    float           scl;
    int             nScl = 10;
    float           scl_range = 1.0;    // (0 - 2.0)
    RecLoop         *scLp, *ch;
    float           lap_w = 0.2;    // 0.3
    BOOL            dbg = YES;

    img = [self copy];
    ch = [img topLoop]; // should check for existence of ch loop
    scLp = [RecLoop loopWithDataLength:nScl];
    imgs = [RecImage imageOfType:RECIMAGE_REAL withLoops:scLp, [img yLoop], [img zLoop], [img xLoop], nil];    // coronal

    for (i = 0; i < nScl; i++) {
        scl = (float)i * scl_range / nScl;
        printf("scale %d = %f\n", i, scl);

        sft_scl = [sft copy];            // pixels
        [sft_scl multByConst:scl];
        pws = [pw correctZShift:sft_scl]; // z shift
        if (i == nScl/2) {
            tmp_img = [pws combineForLoop:ch];
            [tmp_img saveAsKOImage:@"IMG_pws"];
        }
        raw = [pws copy];
    [raw magnitude];    // necessary for akasaka etc
        [raw fft1d:[raw xLoop] direction:REC_INVERSE];
        [raw swapLoop:[raw yLoop] withLoop:[raw zLoop]];
        [grid grid2d:raw to:img];
        img_c = [img combineForLoop:ch];
        [imgs copySlice:img_c atIndex:i forLoop:scLp];
    }
    tmp_img = [imgs copy];
    [tmp_img swapLoop:scLp withLoop:[img yLoop]];
    [tmp_img saveAsKOImage:@"imgs_scl.img"];

// (5) === find best focus scale ====
    printf("5: select best shift\n");
    tmp_img = [imgs copy];
//            [tmp_img laplace2d:REC_FORWARD];
    [tmp_img grad1dForLoop:[tmp_img yLoop]];
    [tmp_img magnitude];
    [tmp_img square];
    [tmp_img gauss3DLP:lap_w];
    [tmp_img swapLoop:scLp withLoop:[img yLoop]];    // ## chk
    [tmp_img saveAsKOImage:@"imgs_scl_lap.img"];

    mxSft = [tmp_img peakIndexForLoop:scLp];    // sharper with larger gradient_sq
    [mxSft saveAsKOImage:@"imgMxSft.img"];        // ok
    [imgs swapLoop:scLp withLoop:[img yLoop]];    // ## chk
    tmp_img = [imgs selectSft:mxSft];            // ### loop ? [scl z x]
    // ####
    [tmp_img saveAsKOImage:@"imgs_final.img"];

    [mxSft gauss2DLP:0.1];
    [mxSft saveAsKOImage:@"imgMxSft_f.img"];
    tmp_img = [imgs selectSftF:mxSft];
    [tmp_img saveAsKOImage:@"imgs_final_F.img"];
    [tmp_img swapLoop:[tmp_img yLoop] withLoop:[tmp_img zLoop]];
    [tmp_img saveAsKOImage:@"imgs_final_ax.img"];
    return imgs;
}

// sft: unit is pixels
// now z only #### make this really 1d
- (RecImage *)correctZShift:(RecImage *)sft;	// z
{
    RecImage    *img, *param;
	float		*p, *q;
	int			i, len;

	// scale sft (pixels -> fracFOV)
	param = [sft copy];
	len = [param xDim];
	[param makeComplex];
	p = [param real];
	q = [param imag];
	// swap real - imag
	for (i = 0; i < len; i++) {
		q[i] = p[i];
		p[i] = 0;
	}
	img = [self ftShiftBy:param]; // ###

    return img;
}

// xy shift (unit ?)
- (RecImage *)correctShift:(RecImage *)sft
{
	return [self ftShiftBy:sft]; // self is not modified
}

// 1D shift for central-view shift estimation
// unit: pixels
- (RecImage *)corr1dWithRef:(RecImage *)ref
{
	RecImage	*sft, *corr;
	RecImage	*tmp1, *ref2;
	float		w = 0.2;	// 0.05
    int         i, mx_i, len;
    float       *p, mx;

	tmp1 = [self copy];
	ref2 = [ref copy];
    // correlation
	corr = [tmp1 xCorrelationWith:ref width:w triFilt:NO];
    [corr saveAsKOImage:@"IMG_corr"];

    tmp1 = [corr avgForLoop:[corr yLoop]];
[tmp1 saveAsKOImage:@"IMG_corr_proj"];
    p = [tmp1 data];
    mx = 0;
    len = [corr xDim];
    for (i = 0; i < len; i++) {
        if (p[i] > mx) {
            mx = p[i];
            mx_i = i;
        }
    }
printf("max at %d (%f)\n", mx_i, mx);
    [corr shift1d:[corr xLoop] by:mx_i - len/2];
//    [corr saveAsKOImage:@"IMG_corr_sft"];

    [corr crop:[corr xLoop] to:20]; // 32
    [corr saveAsKOImage:@"IMG_corr_crop"];   

    // peak detection -> sft
	sft = [corr corrToSft1d:[corr xLoop]];
    mx = [sft meanVal];
    [sft addConst:-mx];
 [sft saveShift:1];
 //[sft saveAsKOImage:@"IMG_sft"];   

	return sft;
}

// move to RecImage.m when done
float
Rec_find_nearest_peak(float *p, int skip, int len)
{
    float       frac;
    int         i, m1, m2, ix, ct = len/2;
    vDSP_Length vix;
    float       mx1, mx2;
    BOOL        incr;

vDSP_maxvi(p, skip, &mx1, &vix, len);   // find max val with index
printf("vDSP: mx = %d %4.3f\n", vix, mx1);

for (i = ct; i < len; i++) {
//    printf("%3.1f ", p[i * skip]);
}
//printf("\n");
for (i = ct; i >= 0; i--) {
//    printf("%3.1f ", p[i * skip]);
}
//printf("\n");
    // find peak nearest to center
    mx1 = mx2 = p[ct * skip];
    m1 = m2 = ct;
    if (p[ct * skip] < p[ct * skip + skip]) {
        incr = YES;
    } else {
        incr = NO;
    }
    for (i = ct+1; i < len; i++) {
        ix = i * skip;
        if (p[ix] > mx1) {
            mx1 = p[ix];
            incr = YES;
        } else {
            if (incr) {
                m1 = i - 1;
                break;  // found m1
            }
            incr = NO;
        }
    }
    if (p[ct * skip] < p[ct * skip - skip]) {
        incr = YES;
    } else {
        incr = NO;
    }
    for (i = ct - 1; i >= 0; i--) {
        ix = i * skip;
        if (p[ix] > mx2) {
            mx2 = p[ix];
            incr = YES;
        } else {
            if (incr) {
                m2 = i + 1;
                break;  // found m1
            }
            incr = NO;
        }
    }
//printf("m1: %d %4.3f,  m2: %d %4.3f\n", m1, mx1, m2, mx2);
    if (mx1 > mx2) {
        ix = m1*skip;
    } else {
        ix = m2*skip;
    }
printf("nearest ix = %d\n", ix);

//mx1 = p[0];
//for (i = 0; i < len; i++) {
//    ix = i * skip;
//    if (p[ix] > mx1) {
//        mx1 = p[ix];
//        m1 = i;
//    }
//}
//printf("sequential mx = %4.3f, at %d\n", mx1, m1);

    if (ix <= 0 || ix >= len * skip) {
        frac = 0;
    } else {
        frac = Rec_find_peak_frac(p + ix, skip, len);
    }
    return ix + frac;
}


- (RecImage *)corrToSft1d:(RecLoop *)lp
{
    void    (^proc)(float *q, float *p, int len, int skip);
 
    proc = ^void(float *q, float *p, int len, int skip) {
    //    *q = Rec_find_peak(p, skip, len);
        *q = Rec_find_nearest_peak(p, skip, len);
    };
	return [self applyProjProc:proc forLoop:lp];
}

// unit: pixels
- (RecImage *)correctShift1d:(RecImage *)sft forLoop:(RecLoop *)lp
{
    RecImage		*img, *param;
	int				len;

	param = [sft copy];
	len = [param xDim];
	img = [self ftShift1d:lp by:param]; // unit: pixels

    return img;
}

// vector delta
- (float)deltaWithPreviousSft:(RecImage *)sft
{
    int     i, n = [self dataLength];
    float   *x, *y, *x0, *y0;
    float   dx, dy, dd, mg;

    if (sft == nil) {   // 1st iteration
        return 1.0;
    }
    x0 = [sft data];
    y0 = x0 + n;
    x = [self data];
    y = x + n;
    dd = mg = 0;
    for (i = 0; i < n; i++) {
        dx = x[i] - x0[i];
        dy = y[i] - y0[i];
        dd += dx*dx + dy*dy;
        mg += x[i]*x[i] + y[i]*y[i];
    }
    dd = sqrt(dd);
    mg = sqrt(mg);

    return dd / mg;
}

// compare magnitude
- (float)deltaWithPreviousSftX:(RecImage *)sft
{
    float   mg, mg0;
    int     i, n = [self dataLength];
    float   *dx, *dy;

    mg0 = 0;
    if (sft != nil) {
        dx = [sft data];
        dy = dx + n;
        for (i = 0; i < n; i++) {
            mg0 += dx[i] * dx[i] + dy[i] * dy[i];
        }
        mg0 = sqrt(mg0);    // calc old mg
    }
    mg = 0; // calc new mg;
    dx = [self data];
    dy = dx + n;
    for (i = 0; i < n; i++) {
        mg += dx[i] * dx[i] + dy[i] * dy[i];
    }
    mg = sqrt(mg);    // calc old mg
    
    return (mg - mg0) / mg;
}

// text dump in pixels
- (void)saveShift:(int)ix
{
	float	*p;
    int     i, j, len;
    int     nPts;
    FILE    *fp;
    char    path[256];

    nPts = [self pixSize]/2;
    sprintf(path, "sft%02d.txt", ix);
    fp = fopen(path, "w");
	p = [self data];
    len = [self xDim];
    for (i = 0; i < len; i++) {
        fprintf(fp, "%d %f\n", i, p[i]);
    }
    fclose(fp);
}

- (void)remove2RR:(RecImage *)pw thres:(float)th	// frac of max
{
	int		i, j, n;
	int		nproj = [pw zDim];
	float	*p, *mgtab, mg, mx, *base;

	[pw fGauss2DLP:0.2];	// window in space
//[pw saveAsKOImage:@"pw_win.img"];

	p = [pw data];
	n = [pw xDim] * [pw yDim];
	mgtab = (float *)malloc(sizeof(float) * nproj);

	for (i = 0; i < nproj; i++) {
		mg = 0;
		for (j = 0; j < n; j++) {
		//	mg += p[j];
			if (mg < p[j]) mg = p[j];
		}
		mg /= n;
		mgtab[i] = mg;
		p += n;
//printf("%f\n", mg);
	}
// baseline removal
	base  = (float *)malloc(sizeof(float) * nproj);
	for (i = 0; i < nproj; i++) {
		base[i] = mgtab[i];
	}
	Rec_smooth(base, nproj, 10);
	for (i = 0; i < nproj; i++) {
	//	printf("%d %f %f\n", i, mgtab[i], mgtab[i] - base[i]);
		mgtab[i] -= base[i];
	}
// normalize
	mx = 0;
	for (i = 0; i < nproj; i++) {
		if (mx < mgtab[i]) mx = mgtab[i];
	}
	for (i = 0; i < nproj; i++) {
		mgtab[i] /= mx;
	}

	p = [self data] + [self dataLength] * 2;	// den
	n = [self xDim];
	for (i = 0; i < nproj; i++) {
		if (mgtab[i] > th) {
			printf("2RR: proj %03d) removed (mg = %4.3f)\n", i, mgtab[i]);
			// set density to 0 (flag)
			for (j = 0; j < n; j++) {
				p[i * n + j] = 0;
			}
		}
	}
	free(mgtab);
	free(base);
}

- (void)removeDeepBreath:(RecImage *)sft thres:(float)th	// frac FOV
{
	int		i, j, n, nproj;
	float	*p, *py, dy;

	nproj = [sft xDim];
	py = [sft data] + [sft dataLength];			// dy

	p = [self data] + [self dataLength] * 2;	// den
	n = [self xDim];
	for (i = 0; i < nproj; i++) {
		dy = py[i];
	//	if (fabs(dy) > th) {
		if (dy > th) {
			printf("DB: proj %03d) removed (dy = %4.3f)\n", i, dy);
			for (j = 0; j < n; j++) {
				p[i * n + j] = 0;
			}
		}
	}
}

// inut is multi-channel (reconstructed) complex image
- (RecImage *)coilMap:(RecLoop *)ch
{
	RecImage	*img_avg, *img_coil;
	RecImage	*mg_a, *mg_c;
	int			i, j, ix, nch, len;
	float		*a, *c, *p, *q, *pp, *qq;
	float		mx;
//	float		mg, ph;

	img_coil = [self copy];
//	[img_coil gauss2DLP:0.05];		// smoothing before div doesn't work ###
//[img_coil saveAsKOImage:@"img_filt.img"];
	img_avg = [img_coil sumForLoop:ch];
[img_avg saveAsKOImage:@"img_avg.img"];

// ### remove small mg pix
	mg_a = [img_avg copy];
	[mg_a magnitude];
	mx = [mg_a maxVal];

	mg_c = [img_coil copy];
	[mg_c magnitude];
	nch = [ch dataLength];
	len = [mg_a dataLength];
	a = [mg_a data];
	c = [mg_c data];

	p = [img_avg data];
	q = p + len;

	pp = [img_coil data];
	qq = pp + [img_coil dataLength];

	for (i = ix = 0; i < nch; i++) {
		for (j = 0; j < len; j++, ix++) {
			if (a[j] < mx * 0.05) {
				pp[ix] = 0.0;
				qq[ix] = 0;
		//	p[j] = mx * 0.01;
		//	q[j] = 0;
		//	if (a[j] != 0) {
		//		mg = c[ix] / a[j];
		//		ph = atan2(pp[ix], qq[ix]) - atan2(p[j], q[j]);
		//		pp[ix] = mg * cos(ph);
			//	qq[ix] = mg * sin(ph);
		//		qq[ix] = 0;
			}
		}
	}

//[img_avg saveAsKOImage:@"img_avg_lmt.img"];

	[img_coil cpxDivImage:img_avg];
	return img_coil;
}

// expand... increase size and fill with outermost value (mostly for interpolation)
- (void)expandLoop:(RecLoop *)lp to:(int)sz
{
	RecLoopControl		*lc;
	int					sz0, lo, hi;
	int					i, j, loopLen, skip, dataLen;
	float				valr, vali, *p, *q;

	sz0 = [lp dataLength];
	lp = [self zeroFill:lp to:sz];
	lc = [self control];
	[lc deactivateLoop:lp];
	loopLen = [lc loopLength];
	skip = [self skipSizeForLoop:lp];
	dataLen = [self dataLength];
	lo = (sz - sz0) / 2;
	hi = sz - lo - 1;
	if ((sz - sz0) % 2 != 0) {
		lo += 1;
	}

	for (i = 0; i < loopLen; i++) {
		p = [self currentDataWithControl:lc];
		q = p + dataLen;
		// lower edge
		valr = p[lo * skip];
		vali = q[lo * skip];
		for (j = 0; j < lo; j++) {
			p[j * skip] = valr;
			q[j * skip] = vali;
		}
		// higher edge
		valr = p[hi * skip];
		vali = q[hi * skip];
		for (j = hi + 1; j < sz; j++) {
			p[j * skip] = valr;
			q[j * skip] = vali;
		}
		[lc increment];
	}
}

- (void)echoFiltForLoop:(RecLoop *)lp	// q&d filter for residual echo removal
{
	[self fGauss1DHP:0.5 forLoop:lp center:-30 frac:0.8];
}

// testing...
// input is raw [ch, z, proj, kx]
- (void)ktPCA
{
	RecImage		*kt;
	int				i, j, tDim, iDim, ix;
	int				xDim, yDim;
	float			*p, *q, *pp, *qq;

	RecImage		*img, *uimg;
	Num_mat			*A, *U, *Vt;
	Num_vec			*s;
	int				nr, nc, n;

//kt = [self combineForLoop:[self topLoop]];
//[kt makeComplex];
	kt = [self copy];
//[kt gauss1DLP:0.2 forLoop:[kt xLoop]];
[kt saveAsKOImage:@"IMG_kt.img"];

//exit(0);
	tDim = [kt yDim];
	iDim = [kt nImages];
	xDim = [kt xDim];
	yDim = [kt yDim];

// === SVD
	//make A matrix
	img = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:iDim yDim:tDim];
	pp = [img real];
	qq = [img imag];
	p = [kt data];
	q = [kt imag];

// ### select k-center line (x == xDim/2)
	for (i = 0; i < tDim; i++) {
		for (j = 0; j < iDim; j++) {
			ix = j * xDim * yDim + i * xDim + xDim/2;
			pp[i * iDim + j] = p[ix];
			qq[i * iDim + j] = q[ix];
		}
	}
	[img saveAsKOImage:@"IMG_A"];

	nr = [img yDim];
	nc = [img xDim];
	n = MIN(nr, nc);
	A = Num_im_to_m(img);
	// make U, s, Vt
	U = Num_new_cmat(nr, n);
	Vt = Num_new_cmat(n, nc);
	s = Num_new_vec(n);
	// call LAPACK
//	Num_svd(A, U, s, Vt);
// ### update svd

	img = Num_v_to_im(s);
	[img saveAsKOImage:@"IMG_s"];

	uimg = Num_m_to_im(U);
	[uimg trans];
	[uimg saveAsKOImage:@"IMG_U"];
	
	img = Num_m_to_im(Vt);
	[img saveAsKOImage:@"IMG_Vt"];
}

- (void)svd_U:(RecImage **)uimg S:(RecImage **)simg Vt:(RecImage **)vtimg
{
	Num_mat			*A, *U, *Vt;
	Num_vec			*s;
	int				nr, nc, n;

	nr = [self yDim];
	nc = [self xDim];
	n = MIN(nr, nc);
	A = Num_im_to_m(self);	// alloc and copy
	// make U, s, Vt
//	if (type == RECIMAGE_COMPLEX) {
	if (pixSize == 2) {	// COMPLEX, MAP etc
		U = Num_new_cmat(nr, n);
		Vt = Num_new_cmat(n, nc);
	} else {
		U = Num_new_mat(nr, n);
		Vt = Num_new_mat(n, nc);
	}
	s = Num_new_vec(n);

	// call LAPACK
//	Num_svd(A, //U, s, Vt);
// ### update svd

	*simg = Num_v_to_im(s);	// alloc and copy
	*uimg = Num_m_to_im(U);
	*vtimg = Num_m_to_im(Vt);

	Num_free_mat(A);
	Num_free_mat(U);
	Num_free_vec(s);
	Num_free_mat(Vt);
}

// self is sft
- (RecImage *)sortZSft
{
	RecImage	*sftTab = [self copy];	// self is not modified
	float		*sft;
	int			len;
//	int			(^compar)(const void *, const void *);

	sft = [sftTab data];
	len = [sftTab dataLength];

	qsort_b_float(sft, len, 1);

	return sftTab;
}

// self is sft (not modified)
// entry is view # in original pw data
- (RecImage *)binTabWithNBins:(int)nBins binSize:(int)sz
{
	RecImage	*binTab = [RecImage imageOfType:RECIMAGE_REAL xDim:sz yDim:nBins];
	RecImage	*sftTab = [self copy];
	int			i, j, ix, st, len = [self xDim];
	float		*p, *q, *sft;
	float		*mn, *mx, *avg;

	p = [self data];
	q = [binTab data];
	sft = [sftTab data];
	mn = (float *)malloc(sizeof(float) * nBins);
	mx = (float *)malloc(sizeof(float) * nBins);
	avg = (float *)malloc(sizeof(float) * nBins);

	sftTab = [self sortZSft];
	sft = [sftTab data];

	// find boundary values
	for (i = 0; i < nBins; i++) {
		st = (len - sz) / (nBins - 1) * i;
		mn[i] = sft[st];
		mx[i] = sft[st + sz - 1];
	//	printf("mn = %f, mx = %f\n", mn[i], mx[i]);
	}

	// make binTab
	for (i = 0; i < nBins; i++) {
		printf("mn/mx = %f/%f\n", mn[i], mx[i]);
		ix = i * sz;
		avg[i] = 0;
		for (j = 0; j < len; j++) {
			if (p[j] >= mn[i] && p[j] <= mx[i]) {
				q[ix++] = j;
				avg[i] += p[j];
			}
		}
		avg[i] /= sz;
		printf("%d %f\n", i, avg[i]);	// -0.5, 0,0.5 (probably wrong)
	}
	free(mn);
	free(mx);
	free(avg);

	return binTab;
}

// block base
// ### sft ok, pw NG ###
/*
- (RecImage *)selectEntryForLoop:(RecLoop *)lp withTab:(RecImage *)tab
{
	RecImage	*img;
	RecLoop		*newLp;
	float		*origIx = [tab data];
	float		*origData = [self data];
	int			origLen = [lp dataLength];
	int			origSkip = [self skipSizeForLoop:lp];
    void		(^proc)(float *p, int len, int skip);


	newLp = [tab xLoop];

	img = [RecImage imageWithImage:self];
	[img replaceLoop:lp withLoop:newLp];

    proc = ^void(float *p, int len, int skip) {
        int     i, ix;
		
		for (i = ix = 0; i < len; i++, ix += skip) {
			p[ix] = origData[(int)(origIx[i]) * origSkip];
		}
    };
//    [img apply1dProc:proc forLoop:newLp]; 

	return img;
}
*/

- (RecImage *)selectEntryForLoopX:(RecLoop *)lp withTab:(RecImage *)tab
{
	RecImage		*img;
	RecLoopControl	*srcLc, *dstLc;
	RecLoopIndex	*px, *npx, *cx;
	RecLoop			*newLp;
	int				prj, nprj, ch, i, k, nProj, nCh, len;
	float			*view;
	float			*p1, *p2;	// src, dst

	// replace pe with xLoop, add yLoop at top level
	newLp = [tab xLoop];
	nProj = [tab xDim];
	view = [tab data];

	srcLc = [self control];
	dstLc = [RecLoopControl controlWithControl:srcLc];
	[dstLc replaceLoop:lp withLoop:newLp];

// make below more general ###

//	pw  [ch, prj, sl, rd]	pwb  [bin, ch, nprj, sl rd]	// #not used
//	raw [ch, sl, prj, rd]	rawb [bin, ch, sl, nprj, rd]
//	trj [prj, rd]			prjb [bin, prj, rd];
//	sft [prj]	// ### not supported => implement more general form

	px = [srcLc loopIndexForLoop:lp];				// proj
	npx = [dstLc loopIndexForLoop:newLp];			// newProj
	if ([srcLc dim] > 3) {	// pw -> raw
		cx = [srcLc loopIndexForLoop:[srcLc topLoop]];
		nCh = [cx dataLength];
		len = [self xDim] * [self yDim];
	} else {
		cx = [RecLoopIndex pointLoopIndex];
		nCh = 1;
		len = [self xDim];
	}
	img = [RecImage imageOfType:[self type] withControl:dstLc];

	// copy data
	for (nprj = 0; nprj < nProj; nprj++) {
		prj = (int)view[nprj]; // view# in src
		[npx setCurrent:nprj];
		[px  setCurrent:prj];
		for (ch = 0; ch < nCh; ch++) {
			[cx setCurrent:ch];
			p1 = [self currentDataWithControl:srcLc];
			p2 = [img currentDataWithControl:dstLc];

			for (k = 0; k < pixSize; k++) {
				for (i = 0; i < len; i++) {
					p2[i] = p1[i];
				}
				p1 += [self dataLength];
				p2 += [img dataLength];
			}
		}
	}

	return img;
}

// raw: copy unit is not xy ###
- (RecImage *)breakLoop:(RecLoop *)lp intoBins:(RecImage *)bTab
{
	RecImage		*img;
	RecLoopControl	*srcLc, *dstLc;
	RecLoopIndex	*bx, *px, *npx, *cx;
	RecLoop			*newLp, *binLp;
	int				bin, prj, nprj, ch, i, k, nBin, nProj, nCh, len;
	float			*view;
	float			*p1, *p2;	// src, dst

	// replace pe with xLoop, add yLoop at top level
	newLp = [bTab xLoop];
	nProj = [bTab xDim];
	binLp = [bTab yLoop];
	nBin = [bTab yDim];
	view = [bTab data];

	srcLc = [self control];
	dstLc = [RecLoopControl controlWithControl:srcLc];
	[dstLc replaceLoop:lp withLoop:newLp];
	[dstLc insertLoop:binLp atIndex:0];

//	pw  [ch, prj, sl, rd]	pwb  [bin, ch, nprj, sl rd]	// #not used
//	raw [ch, sl, prj, rd]	rawb [bin, ch, sl, nprj, rd]
//	trj [prj, rd]			prjb [bin, prj, rd];
	bx = [dstLc loopIndexForLoop:[dstLc topLoop]];	// bin
	px = [srcLc loopIndexForLoop:lp];				// proj
	npx = [dstLc loopIndexForLoop:newLp];			// newProj
	if ([srcLc dim] > 3) {	// pw -> raw
		cx = [srcLc loopIndexForLoop:[srcLc topLoop]];
		nCh = [cx dataLength];
		len = [self xDim] * [self yDim];
	} else {
		cx = [RecLoopIndex pointLoopIndex];
		nCh = 1;
		len = [self xDim];
	}
	img = [RecImage imageOfType:[self type] withControl:dstLc];

	// copy data
	for (bin = 0; bin < nBin; bin++) {
		[bx setCurrent:bin];
		for (nprj = 0; nprj < nProj; nprj++) {
			prj = (int)view[bin * nProj + nprj]; // view# in src
			[npx setCurrent:nprj];
			[px  setCurrent:prj];
		//	printf("bin:prj:nprj = %d:%d:%d\n", bin, prj, nprj); // pw ok, traj ok
			for (ch = 0; ch < nCh; ch++) {
				[cx setCurrent:ch];
				p1 = [self currentDataWithControl:srcLc];
				p2 = [img currentDataWithControl:dstLc];

				for (k = 0; k < pixSize; k++) {
					for (i = 0; i < len; i++) {
						p2[i] = p1[i];
					}
					p1 += [self dataLength];
					p2 += [img dataLength];
				}
			}
		}
	}

	return img;
}

// ### 7-25-2017
// ROI size fixed to 32 x 32 x 32
// covers entire FOV
//  -> make mask and remove background
//
- (RecImage *)defVectorWithRef:(RecImage *)ref
{
	RecImage	*img1, *img2;
	RecImage	*sub1, *sub2, *corr;
	RecImage	*dispTab, *mask, *corrVal;
	float		*dx, *dy, *dz;
	float		*mx, *mn, *cv;
	float		thresC, thresM;
	int			szx, szy, szz;
	int			nx, ny, nz;
	int			xDim, yDim, zDim;
	int			i, j, k, ix, len;
	int			xc, yc, zc;
	RecVector	v;
//	int			sz = 32;

	img1 = [ref copy];
	img2 = [self copy];

	if (0) {
		[img1 subImage:img2];
		[img1 saveAsKOImage:@"diff1.img"];
		exit(0);
	}

// FOV
	xDim = [img1 xDim];
	yDim = [img1 yDim];
	zDim = [img1 zDim];
// ROI size fixed
	szx = 64;
	szy = 64;
	szz = 16;
// # of ROI's
	nx = xDim * 2 / szx - 1;
	ny = yDim * 2 / szy - 1;
	nz = zDim * 2 / szz - 1;

	dispTab = [RecImage imageOfType:RECIMAGE_VECTOR xDim:nx yDim:ny zDim:nz];
	mask = [RecImage imageOfType:RECIMAGE_REAL withImage:dispTab];
	corrVal = [RecImage imageOfType:RECIMAGE_REAL withImage:dispTab];

	len = [dispTab dataLength];
	dx = [dispTab data];
	dy = dx + len;
	dz = dy + len;
	mn = [mask real];
	cv = [corrVal real];

	for (i = ix = 0; i < nz; i++) {
		zc = (i + 1) * szz/2;
		for (j = 0; j < ny; j++) {
			yc = (j + 1) * szy/2;
			for (k = 0; k < nx; k++, ix++) {
				xc = (k + 1) * szx/2;
				sub1 = makeSubImage3(img1, szx, szy, szz, xc, yc, zc);
				sub2 = makeSubImage3(img2, szx, szy, szz, xc, yc, zc);
			//	mn[ix] = [sub1 meanVal];
				mn[ix] = [sub1 rmsVal];
				corr = [sub2 xyzCorrelationWith:sub1 width:0.15 triFilt:YES];
			//	v = [corr findPeak3D];
				v = [corr findPeak3DwithMax:cv + ix];
				// fill displacement tab
				dx[ix] = v.x;
				dy[ix] = v.y;
				dz[ix] = v.z;
			}
		}
	}
[mask saveAsKOImage:@"IMG_ask.img"];
[corrVal saveAsKOImage:@"IMG_cvar.img"];
[dispTab saveAsKOImage:@"IMG_disp.img"];
printf("max corr = %f\n", [corrVal maxVal]);

	thresC = [corrVal maxVal];
	thresC *= 0.6;	// threshold
	thresM = [mask maxVal];
	thresM *= 0.5;
	for (i = 0; i < len; i++) {
	//	if (cv[i] < thresC || mn[i] < thresM) {
		if (mn[i] < thresM) {
			dx[i] = dy[i] = dz[i] = 0;
		}
	}
[dispTab saveAsKOImage:@"dispTab.img"];

	return dispTab;
}

// self: def (3x3x3, 1D, complex), sortTab:[pe]
- (RecImage *)defToSftWithTab:(RecImage *)tab sft:(RecImage *)sft // sft [blk, pe]
{
	RecImage	*bsft, *ssft;
	RecLoop		*bLp, *peLp;
	int			blk, view, nView, ix;
	float		*x, *y, *dx, *dy, *dz, *tabp;
	float		*zsft, scale, th;
	float		z0, z1, a, b;

	peLp = [tab xLoop];
	bLp  = [self xLoop];

	nView = [tab dataLength];

	dx = [self data];
	dy = dx + dataLength;
	dz = dy + dataLength;

	tabp = [tab data];

	// calc zsft amount corresponding to bin-difference (diff)
	// def vector corresponding above diff
	// final shift is def(rotated) x zsft / diff

	// estimate scale (in terms of def vector) of zsft
	ssft = [sft sortZSft];
	zsft = [ssft data];
	z0 = z1 = 0;
	for (view = 0; view < nView/2; view++) {
	//	printf("%d %f\n", view, zsft[view]);
		z0 += zsft[view];
		z1 += zsft[view + nView/2];
	}
	z0 /= nView/2;
	z1 /= nView/2;
	a = 1.0 / (z1 - z0);				// center is mean -> probably end expiratory phase is better
	b = -0.5 * (z1 + z0) / (z1 - z0);
//	printf("z0 / z1 = %f / %f, a = %f, b = %f\n", z0, z1, a, b);

	bsft  = [RecImage imageOfType:RECIMAGE_MAP withLoops:bLp, peLp, nil];
	x = [bsft real];
	y = [bsft imag];
	zsft = [sft data];

	for (blk = 0; blk < 27;	blk++) {
		for (view = 0; view < nView; view++) {
			ix = blk * nView + view;
			scale = zsft[view] * a + b;
		//	printf("%d %f\n", view, scale);
			th = tabp[view];
			// calc sft vector projection to pw plane
			// (rotate def vector in xy plane)
			x[ix] = (dx[blk] * sin(th) + dy[blk] * cos(th)) * scale;
			y[ix] = dz[blk] * scale;
		}
	}

	return bsft;
}

- (RecImage *)outerBlock:(RecBlock3)blk
{
	RecImage		*img;
	RecImage		*kern;	// 3D kernel
	RecLoop			*kx, *ky, *kz;
	float			*p;
	int				i, j, k, ii, jj, kk;
	int				xDim, yDim;

	kx = [self xLoop];
	ky = [self yLoop];
	kz = [self zLoop];

	xDim = [self xDim];
	yDim = [self yDim];

	kern = [RecImage imageOfType:RECIMAGE_REAL withLoops:kz, ky, kx, nil];
	[kern setConst:1.0];
	p = [kern data];
	for (i = 0; i < blk.zW; i++) {
		ii = i + blk.z0;
		for (j = 0; j < blk.yW; j++) {
			jj = j + blk.y0;
			for (k = 0;  k < blk.xW; k++) {
				kk = k + blk.x0;
				p[((ii * yDim) + jj) * xDim + kk] = 0.0;
			}
		}
	}
	img = [self copy];
	[img multByImage:kern];
	return img;
}

- (RecImage *)innerBlock:(RecBlock3)blk
{
	RecImage		*img;
	RecImage		*kern;	// 3D kernel
	RecLoop			*kx, *ky, *kz;
	float			*p;
	int				i, j, k, ii, jj, kk;
	int				xDim, yDim;

	kx = [self xLoop];
	ky = [self yLoop];
	kz = [self zLoop];

	xDim = [self xDim];
	yDim = [self yDim];

	kern = [RecImage imageOfType:RECIMAGE_REAL withLoops:kz, ky, kx, nil];
	p = [kern data];
	for (i = 0; i < blk.zW; i++) {
		ii = i + blk.z0;
		for (j = 0; j < blk.yW; j++) {
			jj = j + blk.y0;
			for (k = 0;  k < blk.xW; k++) {
				kk = k + blk.x0;
				p[((ii * yDim) + jj) * xDim + kk] = 1.0;
			}
		}
	}
	img = [self copy];
	[img multByImage:kern];
	return img;
}

// self is PW
- (RecImage *)projectBlock:(RecBlock3)blk theta:(RecImage *)theta
{
	RecImage	*projImg;
	RecImage	*param;
	RecLoop		*bxLp, *byLp;
	float		x0, y0;
	float		xc, yc;
	int			i;
	int			xDim, yDim, zDim;
	float		*tab;
	float		*sft, th;

	xDim = [self xDim];	// x dim
	yDim = [self yDim];	// z dim
	zDim = [self zDim];	// nProj

	bxLp = [RecLoop loopWithDataLength:blk.xW];
	byLp = [RecLoop loopWithDataLength:blk.zW];

	xc = (blk.x0 + blk.xW/2) - xDim/2;
	yc = (blk.y0 + blk.yW/2) - xDim/2;
	y0 = blk.z0;

	projImg = [self copy];
	param = [RecImage imageOfType:RECIMAGE_MAP withImage:theta];
	sft = [param real];
	tab = [theta data];
	for (i = 0; i < zDim; i++) {
		th = tab[i];
//printf("%d %f\n", i, th);
	//	x0 = - yc * cos(th) - xc * sin(th); // + blk.xW/2;
		x0 =  yc * sin(th) + xc * cos(th); // + blk.xW/2;
		sft[i] = x0;
//printf("%d %f\n", i, x0);
	}
	projImg = [projImg ftShiftBy:param];	// not im-place op
	[projImg replaceLoop:[projImg xLoop] withLoop:bxLp];
	[projImg replaceLoop:[projImg yLoop] withLoop:byLp offset:y0];

	return projImg;
}

- (RecImage *)scaleMapToImage:(RecImage *)img
{
	RecImage	*map = [self copy];
	int			xDim = [img xDim];
	int			yDim = [img yDim];
	int			zDim = [img zDim];
	int			xd = [self xDim];
	int			yd = [self yDim];
	int			zd = [self zDim];

	map = [map scale1dLoop:[map xLoop] by:(float)xDim/(xd + 1) to:xDim];
	map = [map scale1dLoop:[map yLoop] by:(float)yDim/(yd + 1) to:yDim];
	map = [map scale1dLoop:[map zLoop] by:(float)zDim/(zd + 1) to:zDim];

	return map;
}

// updated on 12-27-2019 ... ### not done yet
// input  loops: slc, scale, y, x.
// sftIx  loops: slc,		 y, x.
// output loops: slc,		 y, x.
- (RecImage *)selectSft:(RecImage *)sftIx
{
	RecImage	*img;			// result
	float		*pp;
	RecLoop		*sclLp = [self loopAtIndex:1];
	RecLoop		*slcLp = [self topLoop];
	RecLoop		*xLp = [self xLoop];
	RecLoop		*yLp = [self yLoop];
	float		*p = [self data];
	float		*sft = [sftIx data];
	int			imgSize, nSlc, nScale;
	int			i, k, mix;
	int			ix1, ix2;

	// real only for the moment
	img = [RecImage imageOfType:RECIMAGE_REAL withLoops:slcLp, yLp, xLp, nil];
	pp = [img data];
	imgSize = [img xDim] * [img yDim];
	nSlc = [slcLp dataLength];
	nScale = [sclLp dataLength];
	for (k = 0; k < nSlc; k++) {
		for (i = 0; i < imgSize; i++) {
			ix1 = k * imgSize + i;
			mix = sft[ix1];
			ix2 = ((k * nScale) + mix) * imgSize + i;	// ###
			pp[ix1] = p[ix2];
		}
	}

	return img;
}

// with fraction
- (RecImage *)selectSftF:(RecImage *)sftIx
{
	RecImage	*img;			// result
	float		*pp;
	RecLoop		*sclLp = [self loopAtIndex:1];
	RecLoop		*slcLp = [self topLoop];
	RecLoop		*xLp = [self xLoop];
	RecLoop		*yLp = [self yLoop];
	float		*p = [self data];
	float		*sft = [sftIx data];
	int			imgSize, nSlc, nScale;
	int			i, k, mix;
	float		frac;
	int			ix1, ix2;

	// real only for the moment
	img = [RecImage imageOfType:RECIMAGE_REAL withLoops:slcLp, yLp, xLp, nil];
	pp = [img data];
	imgSize = [img xDim] * [img yDim];
	nSlc = [slcLp dataLength];
	nScale = [sclLp dataLength];
	for (k = 0; k < nSlc; k++) {
		for (i = 0; i < imgSize; i++) {
			ix1 = k * imgSize + i;
			mix = floor(sft[ix1]);
			frac = sft[ix1] - mix;
			ix2 = ((k * nScale) + mix) * imgSize + i;
			pp[ix1] = p[ix2] * (1 - frac) + p[ix2 + imgSize] * frac;
		}
	}

	return img;
}

// self is sft
// in current recon scheme, scale of sft does not have significanse
- (void)scaleToBins
{
	RecImage	*sft;	// sorted copy of self
	int			i, n, n2;
	float		*p;
	float		mn0, mn1;

	sft = [self sortZSft];
	p = [sft data];
	n = [sft dataLength];
	mn0 = mn1 = 0;
	n2 = n/2;
	for (i = 0; i < n2; i++) {
		mn0 += p[i + n2];
		mn1 += p[i];
	}
	mn0 /= n2;
	mn1 /= n2;
	printf("=== mn0 = %f, mn1 = %f\n", mn0, mn1);

	p = [self data];
	for (i = 0; i < n; i++) {
		p[i] = (p[i] - mn0) / (mn1 - mn0);
	//	printf("%d %f\n", i, p[i]);
	}
}

// output is self
- (void)scaleVector:(RecVector)v shift:(RecImage *)sft rotate:(RecImage *)tab
{
	int		i, n = [self dataLength];	// view
	float	*th = [tab real];
	float	*z = [sft real];	// b0-b1 distance
	float	*x = [self real];	// pixels (out)
	float	*y = [self imag];	// pixels (out)

	for (i = 0; i < n; i++) {
		y[i] = z[i] * v.z * z[i];
		x[i] = z[i] * (v.x * cos(th[i]) - v.y * sin(th[i]));
	}
}

// #####
// *reproduce s-o-s
// *add mag calc / thres
// *div by weight
- (RecImage *)combineForLoop:(RecLoop *)lp threshold:(float)th
{
    RecImage    *img;
    void    (^proc)(float *qr, float *qi, float *pr, float *pi, int llen, int lskip, int rlen, int rskip);
 
    proc = ^void(float *qr, float *qi, float *pr, float *pi, int llen, int lskip, int rlen, int rskip) {
        int     i, j, ix;
        float	sumr, sumi;
		float	mg, mx, wt;
        for (i = 0; i < rlen; i++) {
			mx = 0;
			for (j = ix = 0; j < llen; j++, ix += lskip) {
				mg = sqrt(pr[ix] * pr[ix] + pi[ix] * pi[ix]);
				if (mx < mg) mx = mg;
			}
			mx *= th;
		//	mx *= 0.1;
			
			sumr = sumi = wt = 0;
			for (j = ix = 0; j < llen; j++, ix += lskip) {
				mg = sqrt(pr[ix] * pr[ix] + pi[ix] * pi[ix]);
				if (mg > mx) {
					sumr += pr[ix] * pr[ix];
					sumi += pi[ix] * pi[ix];
					wt += mg;
				}
			}
        //    *qr = sqrt(sumr + sumi);
            *qr = (sumr + sumi) / wt;
            qr += rskip;
            pr += rskip;
            pi += rskip;
        }
    };
	img = [self applyCombineProc:proc forLoop:lp];
    [img takeRealPart];
    return img;
}

// ===== dbg =====

- (void)saveComb:(RecLoop *)ch as:(NSString *)path
{
	RecImage	*img;

	img = [self copy];
	if (ch) {
		img = [img combineForLoop:ch];
	}
	[img saveAsKOImage:path];
}

- (void)printDefAt:(RecVector)v
{
	int		x, y, z, ix;
	float	*p = [self data];
	int		len = [self dataLength];

	x = v.x;
	y = v.y;
	z = v.z;

	ix = (z * [self yDim] + y) * [self xDim] + x;
	printf("sft = %f/%f/%f\n", p[ix], p[ix + len], p[ix + len*2]);
}

- (void)setMarker:(float)maxVal
{
	[self data][0] = maxVal;
}

- (void)dumpShift
{
	int		i, n = [self xDim];
	float	*p, *q;
	float	*srt, *rng;

	if ([self type] == RECIMAGE_REAL) {
		// sft
		p = [self real];
		// accum
		srt = (float *)malloc(sizeof(float) * n);
		for (i = 0; i < n; i++) {
			srt[i] = p[i];
		}
		qsort_b_float(srt, n, 1);
		
		// range
		rng = (float *)malloc(sizeof(float) * n);
		for (i = 0; i < n/2; i++) {
			rng[i] = srt[i + n/2] - srt[i];
		}
		for (; i < n; i++) {
			rng[i] = 0;
		}

		// plot
		for (i = 0; i < n; i++) {
			printf("%d %f %f %f\n", i, p[i], srt[i], rng[i]);
		}
		free(srt);
		free(rng);
	} else {
		p = [self real];
		q = [self imag];
		for (i = 0; i < n; i++) {
			printf("%d %f %f\n", i, p[i], q[i]);
		}
	}
}

// move to RecKit when done
- (void)normalize2ForLoop:(RecLoop *)lp withMask:(RecImage *)msk	// L2 norm
{
	RecLoopControl	*lc;
	float			*p, *q, mg;
	float			*mk;
	int				i, k, ix;
	int				skip, len, loopLen;

	lc = [self control];
	[lc deactivateLoop:lp];
	loopLen = [lc loopLength];
	len = [lp dataLength];
	skip = [self skipSizeForLoop:lp];
	mk = [msk data];

	if ([self type] == RECIMAGE_COMPLEX) {
		for (k = 0; k < loopLen; k++) {
			mg = 0;
			p = [self currentDataWithControl:lc];
			q = p + [self dataLength];
			for (i = ix = 0; i < len; i++, ix += skip) {
				mg += p[ix]*p[ix] + q[ix]*q[ix];
			}
			if (mg != 0) {
				mg = mk[k] / sqrt(mg);
				for (i = ix = 0; i < len; i++, ix += skip) {
					p[ix] *= mg;
					q[ix] *= mg;
				}
			}
			[lc increment];
		}
	} else {
		for (k = 0; k < loopLen; k++) {
			mg = 0;
			p = [self currentDataWithControl:lc];
			for (i = ix = 0; i < len; i++, ix += skip) {
				mg += p[ix]*p[ix];
			}
			if (mg != 0) {
				mg =mk[k] / sqrt(mg);
				for (i = ix = 0; i < len; i++, ix += skip) {
					p[ix] *= mg;
				}
			}
			[lc increment];
		}
	}
}

// ==== testing... ======
// mag filter
- (void)magLPF:(float)w
{
	RecImage	*mg, *phs;

	mg = [self copy];
	[mg magnitude];
	phs = [self copy];
	[phs phase];
	[mg gauss2DLP:w];
	[self copyImage:mg];
	[self makeComplexWithPhs:phs];
}

// polar filter
- (void)polarLPF:(float)w
{
	RecImage	*mg, *phs;

	mg = [self copy];
	[mg magnitude];
	[mg gauss2DLP:w];
	phs = [self copy];
//	[phs phase];
	[phs toUnitImage];
	[phs gauss2DLP:w];
	[phs atan2];

	[self copyImage:mg];
	[self makeComplexWithPhs:phs];
}

//===========

@end

// ===== c func =====

void
Rec_calcTransformParam(float *v, float *x, float *y, float *dx, float *dy, int nfree)
{
//    float           *mat, *mp, *aplus, *u;
    Num_mat         *mat, *aplus;
    Num_vec         *uvec, *vvec;
    float           *mp;
    float           sumx, sumy;
    int             i;

    switch (nfree) {
    case 1 :    // sft
        sumx = sumy = 0;
        for (i = 0; i < 4; i++) {
            sumx += dx[i];
            sumy += dy[i];
        }
        v[0] = sumx / 4;
        v[1] = sumy / 4;
        break;
    case 2 :    // sft + scale
    // least suqares solution
        mat     = Num_new_mat(8, 4);
        aplus   = Num_new_mat(4, 8);
        uvec    = Num_new_vec(8);
        vvec    = Num_new_vec(4);
    
        mp = mat->data;
        for (i = 0; i < 4; i++) {
            mp[0] = x[i];
            mp[1] = 0;
            mp[2] = 1;
            mp[3] = 0;
            mp += 4;
        }
        for (i = 0; i < 4; i++) {
            mp[0] = 0;
            mp[1] = y[i];
            mp[2] = 0;
            mp[3] = 1;
            mp += 4;
        }
        // setup u
        for (i = 0; i < 4; i++) {
            uvec->data[i]     = x[i] + dx[i];
            uvec->data[i + 4] = y[i] + dy[i];
        }
        // param is inv(m) * v
//        Num_pinv(aplus, mat);	// ### update using Class version
        Num_mvmul(vvec, aplus, uvec);
        for (i = 0; i < 4; i++) {
            v[i] = vvec->data[i];
        }
        Num_free_mat(mat);
        Num_free_mat(aplus);
        Num_free_vec(uvec);
        Num_free_vec(vvec);
        break;
    case 3 :    // affine
        // least squares solution using all 4 points
        mat     = Num_new_mat(8, 6);
        aplus   = Num_new_mat(6, 8);
        uvec    = Num_new_vec(8);
        vvec    = Num_new_vec(6);
        // setup m
        mp = mat->data;
        for (i = 0; i < 4; i++) {
            mp[0] = x[i];
            mp[1] = y[i];
            mp[2] = 1.0;
            mp[3] = mp[4] = mp[5] = 0.0;
            mp += 6;
        }
        for (i = 0; i < 4; i++) {
            mp[0] = mp[1] = mp[2] = 0.0;
            mp[3] = x[i];
            mp[4] = y[i];
            mp[5] = 1.0;
            mp += 6;
        }
        // setup v
        for (i = 0; i < 4; i++) {
            uvec->data[i]     = x[i] + dx[i];
            uvec->data[i + 4] = y[i] + dy[i];
        }
        // param is inv(m) * v
//        Num_pinv(aplus, mat);	// ### update with Class version
        Num_mvmul(vvec, aplus, uvec);
        for (i = 0; i < 6; i++) {
            v[i] = vvec->data[i];
        }
        Num_free_mat(mat);
        Num_free_mat(aplus);
        Num_free_vec(uvec);
        Num_free_vec(vvec);
        break;
    case 4 :    // homography
        mat     = Num_new_mat(8, 8);
        aplus   = Num_new_mat(8, 8);
        uvec    = Num_new_vec(8);
        vvec    = Num_new_vec(8);
        // setup m
        mp = mat->data;
        for (i = 0; i < 4; i++) {
            mp[0] = x[i];
            mp[1] = y[i];
            mp[2] = 1.0;
            mp[3] = mp[4] = mp[5] = 0.0;
            mp[6] = - x[i] * (x[i] + dx[i]);
            mp[7] = - y[i] * (x[i] + dx[i]);
            mp += 8;
        }
        for (i = 0; i < 4; i++) {
            mp[0] = mp[1] = mp[2] = 0.0;
            mp[3] = x[i];
            mp[4] = y[i];
            mp[5] = 1.0;
            mp[6] = - x[i] * (y[i] + dy[i]);
            mp[7] = - y[i] * (y[i] + dy[i]);
            mp += 8;
        }
        // setup v
        for (i = 0; i < 4; i++) {
            uvec->data[i]     = x[i] + dx[i];
            uvec->data[i + 4] = y[i] + dy[i];
        }
        // param is inv(m) * v
        Num_gaussj(mat, NULL);
        Num_mvmul(vvec, aplus, uvec);
        for (i = 0; i < 8; i++) {
            v[i] = vvec->data[i];
        }
        Num_free_mat(mat);
        Num_free_vec(uvec);
        Num_free_vec(vvec);
        break;
    }
}

int
takeROI2d(RecImage *src, RecImage *dst, int xc, int yc, int slice)
{
	int		srcXDim = [src xDim];
	int		srcYDim = [src yDim];
	int		dstXDim = [dst xDim];
	int		dstYDim = [dst yDim];
    float   *srcp = [src data] + slice * srcXDim * srcYDim;
    float   *dstp = [dst data];
    int     i, j;
    int     ii, jj;

    [dst clear];
	for (ii = 0; ii < dstYDim; ii++) {
		i = yc + ii - dstYDim/2;
		if (i < 0 || i >= srcYDim) continue;
		for (jj = 0; jj < dstXDim; jj++) {
			j = xc + jj - dstXDim/2;
			if (j < 0 || j >= srcXDim) continue;
			dstp[ii * dstXDim + jj] = srcp[i * srcXDim + j];
		}
	}
    return 0;
}

RecImage *
makeSubImage2(RecImage *img, int xsz, int ysz, int xc, int yc)
{
    Rec2DRef    *ref;
    int         x0, y0;

    x0 = xc - xsz/2;
    y0 = yc - ysz/2;
    ref = [Rec2DRef refForImage:img];
    [ref setX:x0 y:y0 nX:xsz nY:ysz];
    return [ref makeImage];
}

RecImage *
makeSubImage3(RecImage *img, int xsz, int ysz, int zsz, int xc, int yc, int zc)
{
    Rec3DRef    *ref;
    int         x0, y0, z0;

    x0 = xc - xsz/2;
    y0 = yc - ysz/2;
	z0 = zc - zsz/2;
    ref = [Rec3DRef refForImage:img];
    [ref setX:x0 y:y0 z:z0 nX:xsz nY:ysz nZ:zsz];
    return [ref makeImage];
}

RecImage *
dispToMap2d(RecImage *disp)
{
    RecImage    *map;
    int         i, j, k, ix;
    float       x, y;
    int         len = [disp dataLength];
    int         xDim = [disp xDim];
    int         yDim = [disp yDim];
    int         zDim = [disp zDim];
    float       *dx, *dy;
    float       *mx, *my;

    dx = [disp data];
    dy = dx + len;
    map = [RecImage imageWithImage:disp];
    mx = [map data];
    my = mx + len;
    ix = 0;
    for (k = 0; k < zDim; k++) {
        for (i = 0; i < yDim; i++) {
            y = ((float)i - yDim/2) / yDim;
            for (j = 0; j < xDim; j++, ix++) {
                x = ((float)j - xDim/2) / xDim;
                mx[ix] = dx[ix] / xDim + x;
                my[ix] = dy[ix] / yDim + y;
            }
        }
    }
    return map;
}

// same as dispToMap in lungDef
RecImage *
dispToMap3d(RecImage *disp)
{
    RecImage    *map;
    int         i, j, k, ix;
    float       x, y, z;
    int         len = [disp dataLength];
    int         xDim = [disp xDim];
    int         yDim = [disp yDim];
    int         zDim = [disp zDim];
    float       *dx, *dy, *dz;
    float       *mx, *my, *mz;

    dx = [disp data];
    dy = dx + len;
    dz = dy + len;
    map = [RecImage imageWithImage:disp];
    mx = [map data];
    my = mx + len;
    mz = my + len;
    ix = 0;
    for (k = 0; k < zDim; k++) {
        z = ((float)k - zDim/2) / zDim;
        for (i = 0; i < yDim; i++) {
            y = ((float)i - yDim/2) / yDim;
            for (j = 0; j < xDim; j++, ix++) {
                x = ((float)j - xDim/2) / xDim;
                mx[ix] = -dx[ix] / xDim + x;
                my[ix] = -dy[ix] / yDim + y;
                mz[ix] = -dz[ix] / zDim + z;
            }
        }
    }
    return map;
}

// 1-d phase map
// x:proj, y:slice
// input (pw) is paddle-wheel image (ch, proj, slc, read)
RecImage *
Rec_radDelayMap(RecImage *pw)
{
	float			*re, *im, *p, phs;
	int				i, j, n, len, ix;
	int				xDim, yDim;
	RecImage		*map;
	RecLoopControl	*lc;
	RecLoop			*ch;

	ch = [RecLoop findLoop:@"Channel"];
	if ([pw containsLoop:ch] && [ch dataLength] > 1) {
		pw = [pw avgForLoop:ch];
	}

	xDim = [pw zDim];
	yDim = [pw yDim];
	n = [pw xDim];
	len = [pw dataLength];
	lc = [pw control];
	map = [RecImage imageOfType:RECIMAGE_REAL xDim:xDim yDim:yDim];
	p = [map data];

	[lc rewind];
	[lc deactivateX];
	for (j = 0; j < xDim; j++) {
		for (i = 0; i < yDim; i++) {
			re = [pw currentDataWithControl:lc];
			im = re + len;
			phs = Rec_est_1st(re, im, 1, n);
			ix = i * xDim + j;
			p[ix] = phs * 1000;
			[lc increment];
		}
	}

	return map;
}

RecImage *
Rec_goldenAngleTab(RecLoop *pe, int nRef)
{
	RecImage		*tab;
	int				i, nProj;
	float			*th;

	tab = [RecImage imageOfType:RECIMAGE_REAL withLoops:pe, nil];
	nProj = [tab dataLength];
	th = [tab data];

    for (i = 0; i < nProj; i++) {
		if (i < nProj - nRef) {
			th[i] = Rec_golden_angle(i);
		} else { // phase correction view
			th[i] = (i - nProj + nRef) * M_PI * 2 / nRef;
		}
    }
	return tab;
}

RecImage *
Rec_goldenAngleSort(RecImage *gaTab)
{
	RecImage	*tab;
	float		*src, *dst, *dstIx;
	int			i, len = [gaTab dataLength];
	gaEnt		*sortTab;
	int			(^compar)(const void *, const void *);

	compar = ^int(const void *p1, const void *p2) {
		float pp = ((gaEnt *)p1)->th;
		float qq = ((gaEnt *)p2)->th;
		if (pp > qq) {
			return 1;
		} else 
		if (pp == qq) {
			return 0;
		} else {
			return -1;
		}
	};
	sortTab = (gaEnt *)malloc(sizeof(gaEnt) * len);
	src = [gaTab data];
	for (i = 0; i < len; i++) {
		sortTab[i].th = src[i];
		sortTab[i].origIx = i;
	}
	qsort_b(sortTab, len, sizeof(gaEnt), compar);

	tab = [RecImage imageOfType:RECIMAGE_COMPLEX withImage:gaTab];
	dst = [tab real];
	dstIx = [tab imag];
	for (i = 0; i < len; i++) {
		dst[i]   = sortTab[i].th;
		dstIx[i] = sortTab[i].origIx;
	//	printf("%f %f\n", dst[i], dstIx[i]);
	}

	return tab;
}

void
qsort_b_float(float *p, int n, int order)	// order 1:acsending, -1:descending
{
	int			(^compar)(const void *, const void *);

	// sort sftTab
	compar = ^int(const void *p1, const void *p2) {
		float pp = ((float *)p1)[0];
		float qq = ((float *)p2)[0];
		if (order > 0) {
			if (pp > qq) {
				return 1;
			} else 
			if (pp == qq) {
				return 0;
			} else {
				return -1;
			}
		} else {
			if (pp < qq) {
				return 1;
			} else 
			if (pp == qq) {
				return 0;
			} else {
				return -1;
			}
		}
	};
	qsort_b(p, n, sizeof(float), compar);
}

// not necessary if block sort tab is present in header
void
block_sort(RecImage *raw)
{
	int				i, j, n;
	int				ix, skip;
	int				*tab;
	RecLoop			*lp;
	RecLoopControl	*lc;
	int				loopLen;
    int             dataLength;
	float			*buf, *p;

    n = [raw zDim];
    lp = [raw zLoop];
    dataLength = [raw dataLength];
	buf = (float *)malloc(sizeof(float) * n);
	tab = (int *)malloc(sizeof(int) * n);

    for (i = 0; i < n; i++) {
        if (i == 0) {
            ix = 0;
        } else
        if (i <= n/2) {
            ix = n - i * 2 + 1;
        } else {
            ix = (i - n/2) * 2;
        }
        tab[i] = ix;
    //    printf("%d,", ix);
    }
//    printf("\n");

//    for (i = 0; i < n; i++) {
//        printf("%d\n", tab[i]);
//    }

	skip = [raw skipSizeForLoop:lp];
	lc = [raw control];
	[lc deactivateLoop:lp];
	loopLen = [lc loopLength];
	for (i = 0; i < loopLen; i++) {
		// real
		p = [raw currentDataWithControl:lc];
		for (j = ix = 0; j < n; j++, ix += skip) {
		//	buf[tab[j]] = p[ix];
			buf[j] = p[ix];
		}
		for (j = ix = 0; j < n; j++, ix += skip) {
		//	p[ix] = buf[j];
			p[ix] = buf[tab[j]];
		}
		// imag
		p += dataLength;
		for (j = ix = 0; j < n; j++, ix += skip) {
		//	buf[tab[j]] = p[ix];
			buf[j] = p[ix];
		}
		for (j = ix = 0; j < n; j++, ix += skip) {
		//	p[ix] = buf[j];
			p[ix] = buf[tab[j]];
		}
		[lc increment];
	}
	free(buf);
    free(tab);
}

// amount of z-shift (relative), unit is mm
RecImage *
Rec_load_edge(NSString *path)
{
	RecImage		*sft;
	FILE			*fp;
	const char		*cPath = [path UTF8String];
	char			*buf = (char *)malloc(256);
	char			*sts;
	int				i, nview;
	int				isft;
	float			*sftp, mn;

// open file
	fp = fopen(cPath, "r");
	if (fp == NULL) {
		printf("edge file not found\n");
		exit(0);
	}

	for (i = 0; ; i++) {
		sts = fgets(buf, 256, fp);
		if (sts == NULL) break;
	}
	nview = i;
//	sft = (float *)malloc(sizeof(float) * nview);
	sft = [RecImage imageOfType:RECIMAGE_REAL xDim:nview];
	sftp = [sft data];

	rewind(fp);
	for (i = 0; i < nview; i++) {
        sts = fgets(buf, 256, fp);
		if (sts == NULL) break;
		sscanf(buf, "[new Thread] - correctEdgePoint..%d", &isft);
		sftp[i] = isft * 0.46875;
	}
	fclose(fp);
	free(buf);

	mn = 0;
	for (i = 0; i < nview; i++) {
		mn += sftp[i];
	}
	mn /= nview;
	for (i = 0; i < nview; i++) {
		sftp[i] -= mn;
	}

	return sft;
}

