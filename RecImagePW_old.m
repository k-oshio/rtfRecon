

// ==== OLD =======
// === not used any more ===###
// in/out real
//
// assumption: cross-corr(total) = auto-corr(static) + cross-corr(moving)
//      -> cross-corr(moving) = cross-corr(total) - auto-corr(static)
//

// 4 referelnce locations
// pixels shift from center
int nLoc = 4;
NSPoint hrtCtr = {27, -9}; //{27, -15};
NSRect  pts[4] = {
    {0, -20, 20, 10},   //  0 U
    {0,  20, 20, 10},   //  1 D
    {-30, 0, 40, 30},   //  2 L
    { 30, 0, 40, 30}    //  3 R
};
int w = 5; // 10

- (NSPoint)centerLocForPt:(int)ix view:(int)pe nAng:(int)nView
{
    NSPoint pt;
    float   th, cs, sn;

    th = pe * M_PI / nView;
    cs = cos(th);
    sn = sin(th);

    pt = pts[ix].origin;
    pt.x += [self xDim] / 2 + hrtCtr.x * cs + hrtCtr.y * sn;
    pt.y += [self yDim] / 2;

    return pt;
}

//
// ### not working yet... update range check
// (range check should be there, but corr outside of image boundary is meaningless)
// (check at calling side, too)
//
- (RecImage *)detectShift4WithRef:(RecImage *)temp
{
    RecImage        *corr;      // correlation
    RecImage        *sft;       // shift x 4pts
    NSMutableArray  *loopArray;
    RecLoopControl  *lc1, *lc2;
    RecLoop         *xLp, *yLp, *zLp;
    int             width = w * 2 + 1;
    int             k, nPe, xDim, yDim;
    int             loc, sftLength;
    float           *px, *py;
    Rec3DRef        *ref1, *ref2;
    NSRect          r;
    NSPoint         pos, ctr;
    int             x0, y0; // upper-left corner
    BOOL            mark = NO;

// edge + rho + sideLobe
    [self pwFilter];
    [temp pwFilter];

//##########
// copy, mark corr area, and save for checking
// leave original unmarked (for processing)
    loopArray = [NSMutableArray arrayWithArray:[self loops]];
    [loopArray removeLastObject];
    [loopArray removeLastObject];
    sft = [RecImage imageOfType:RECIMAGE_HOMOG withLoopArray:loopArray];
    sftLength = [sft dataLength];
    xLp = [RecLoop loopWithDataLength:width];
    yLp = [RecLoop loopWithDataLength:width];
    zLp = [RecLoop pointLoop];
    corr = [RecImage imageOfType:RECIMAGE_REAL withLoops:zLp, yLp, xLp, nil];

// proj loop
    lc1 = [RecLoopControl controlForImage:self];
    lc2 = [RecLoopControl controlWithControl:lc1 forImage:corr];
    [lc1 deactivateXY];
    nPe = [lc1 loopLength];
    xDim = [self xDim];
    yDim = [self yDim];

    for (k = 0; k < nPe; k++) {
        for (loc = 0; loc < 4; loc++) {
            // set loc
            ctr = [self centerLocForPt:loc view:k nAng:nPe];
            r.origin.x = ctr.x - pts[loc].size.width/2;
            r.origin.y = ctr.y - pts[loc].size.height/2;
            r.size = pts[loc].size;
            x0 = r.origin.x;
            y0 = r.origin.y;

            ref1 = [Rec3DRef refForImage:self control:lc1]; // target
            ref2 = [Rec3DRef refForImage:temp control:lc1]; // template (fixed)
            [ref1 setX:r.origin.x y:r.origin.y z:0
                nX:pts[loc].size.width nY:pts[loc].size.height nZ:1];
            [ref2 setX:r.origin.x y:r.origin.y z:0
                nX:pts[loc].size.width nY:pts[loc].size.height nZ:1];
            if (mark) {
                [ref1 markWith:-0.8];
                [ref2 markWith:-0.8];
            } else {
                [ref1 normalizedCorrelationWith:ref2 result:corr];
                pos = [corr findPeak2D];
            }
            px = [sft data] + sftLength * loc * 2;
            py = px + sftLength;
            px[k] = - (pos.x - width/2);
            py[k] = - (pos.y - width/2);
        }
        [lc1 increment];
    }
    [self saveAsKOImage:@"pw_tgt.img"];
    [temp saveAsKOImage:@"pw_tmp.img"];
//exit(0);

    return sft;
}

// calc corr for 4 locations, and return sft amount as HOMOG file (float x 8)
- (RecImage *)detectShift4X:(RecImage *)temp
{
    RecImage        *corr;      // correlation
    RecImage        *sft;       // shift x 4pts
    RecLoopControl  *lc1, *lc2;
    RecLoop         *xLp, *yLp;
    NSMutableArray  *loopArray;
    int             width = w*2 + 1;
    int             i, j, k, nPe, xDim, yDim;
    int             loc, sftLength;
    float           *p, *px, *py;
    Rec2DRef        *ref1, *ref2;
    NSRect          r;
    NSPoint         pos, ctr;
    int             x0, y0; // upper-left corner

// edge + rho + sideLobe
    [self pwFilter];
    [temp pwFilter];
    [self saveAsKOImage:@"pw_tgt.img"];
    [temp saveAsKOImage:@"pw_tmp.img"];

    loopArray = [NSMutableArray arrayWithArray:[self loops]];
    [loopArray removeLastObject];
    [loopArray removeLastObject];
    sft = [RecImage imageOfType:RECIMAGE_HOMOG withLoopArray:loopArray];
    sftLength = [sft dataLength];
    xLp = [RecLoop loopWithDataLength:width];
    yLp = [RecLoop loopWithDataLength:width];
    [loopArray addObject:yLp];
    [loopArray addObject:xLp];
    corr = [RecImage imageOfType:RECIMAGE_REAL withLoopArray:loopArray];

// proj loop
    lc1 = [RecLoopControl controlForImage:self];
    lc2 = [RecLoopControl controlWithControl:lc1 forImage:corr];
    [lc1 deactivateXY];
    nPe = [lc1 loopLength];
    xDim = [self xDim];
    yDim = [self yDim];

    for (k = 0; k < nPe; k++) {
        for (loc = 0; loc < 4; loc++) {
            // set loc
            ctr = [self centerLocForPt:loc view:k nAng:nPe];
            r.origin.x = ctr.x - pts[loc].size.width/2;
            r.origin.y = ctr.y - pts[loc].size.height/2;
            r.size = pts[loc].size;
            x0 = r.origin.x;
            y0 = r.origin.y;

            ref1 = [Rec2DRef refForImage:self control:lc1]; // target
            ref2 = [Rec2DRef refForImage:temp control:lc1]; // template (fixed)
            [ref2 setSelection:r];

            p = [corr currentDataWithControl:lc2];
            for (i = 0; i < width; i++) {
                r.origin.y = y0 - w + i;
                for (j = 0; j < width; j++) {
                    r.origin.x = x0 - w + j;
                    [ref1 setSelection:r];
                    p[i * width + j] = [ref1 normalizedCorrelationWith:ref2];
                }
            }
            // calc sft
            pos = Rec_find_peak2(p, width, width);
            px = [sft data] + sftLength * loc * 2;
            py = px + sftLength;
            px[k] = pos.x - width/2;
            py[k] = pos.y - width/2;
        }
        [lc1 increment];
    }
    return sft;
}

- (RecImage *)correctShift4:(RecImage *)sft mark:(BOOL)mk over:(int)ovr
{
    RecImage        *img;    // result
    RecImage        *img_e;
    RecImage        *map;
    RecImage        *param;
    RecLoop         *lp;    // pe
    RecLoopControl  *lc, *pLc;
    int             i, j, n, len, xDim, yDim;
    float           x[4], y[4];
    float           xsft[4], ysft[4];
    float           *v;
    float           *s, *pr;    // data ptr for image, sft, param
    int             nfree = 3;      // 1: sft, 2:sft + scale, 3:affine, 4:homography
    NSPoint         ctr;

    lp = [sft topLoop]; // should be OK for multi-coil case, too
    lc = [RecLoopControl controlForImage:self];
    [lc deactivateXY];
    n = [lp dataLength];
    xDim = [self xDim];
    yDim = [self yDim];

    param = [RecImage imageOfType:RECIMAGE_HOMOG withLoops:lp, nil];
    pLc = [RecLoopControl controlWithControl:lc forImage:param];

//    m = (float *)malloc(8 * 8 * sizeof(float));
    v = (float *)malloc(8 * sizeof(float));

    len = [sft dataLength];

// extract sft params (sft data is always single-channel)
    for (i = 0; i < n; i++) {   // for projections
        s = [sft data];
        for (j = 0; j < nLoc; j++) {    // for 4 locs
            // (x, y) is original locacion
            // (x + xsft, y + ysft) is new location
            xsft[j] = s[i];
            s += len;
            ysft[j] = s[i];
            s += len;
            ctr = [self centerLocForPt:j view:i nAng:n];
            x[j] = ctr.x;
            y[j] = ctr.y;
        }
        // first set mark, then calc param for [-0.5..0.5] to [-0.5..0.5]
        // scale x/y/xsft/ysft to [-0.5 .. 0.5]
        for (j = 0; j < 4; j++) {
            x[j] = (x[j] - xDim/2) / xDim;
            y[j] = (y[j] - yDim/2) / yDim;
            xsft[j] /= xDim;
            ysft[j] /= yDim;
        }
    // ===== make func for 1:shift, 2:shift+scale, 3:affine, 4:homography
        Rec_calcTransformParam(v, x, y, xsft, ysft, nfree);
        // copy v to param
        pr = [param currentDataWithControl:pLc];
        for (j = 0; j < nfree * 2; j++) {
            *pr = v[j];
            pr += n;
        }
        [lc increment];
    }

    if (mk) {
        [self markPoints:sft];
    }
    if (ovr > 0) {  // not correct ...  map has already been calculated with old self !!!
        // wrap-pad
        img_e = [self expandSliceBy:ovr];
    } else {
        img_e = self;
    }
    img = [RecImage imageWithImage:img_e];

// correction
    switch (nfree) {
    case 1 :    // shift
        map = [img mapForShift:param];
        break;
    case 2 :
        map = [img mapForShiftScale:param];
        break;
    case 3 :
        map = [img mapForAffine:param];
        break;
    case 4 :    // homography
        map = [img mapForHomog:param];
        break;
    }
    [img resample:img_e withMap:map];
    if (ovr > 0) {  // not correct ...  map has already been calculated with old self !!!
        [img cropBy:ovr];
    }

    free(v);

    return img;
}

// sft: [pe]
// self:[(ch), pe, sl, rd]
- (void)markPoints:(RecImage *)sft
{
    float           xsft[4], ysft[4];
    int             xDim, yDim;
    int             i, j, n, len;
    int             xi, yi;
    RecLoopControl  *lc, *sLc;
    float           *p, *s;
    float           mark = 30.0;
    RecImage        *pw_c;
    RecLoop         *ch, *pe;
    RecLoopIndex    *peIx;
    NSPoint         ctr;
    int             nView;

// mark all channels
    ch = [RecLoop findLoop:@"Channel"];
    pe = [RecLoop findLoop:@"Phase"];
    lc  = [RecLoopControl controlForImage:self];
    sLc = [RecLoopControl controlWithControl:lc forImage:sft];
    [lc deactivateXY];
    peIx = [lc loopIndexForLoop:pe];

    xDim = [self xDim];
    yDim = [self yDim];
    nView = [pe dataLength];

    n = [lc loopLength];
    len = [sft dataLength];
    for (i = 0; i < n; i++) {
        p = [self currentDataWithControl:lc];
        s = [sft currentDataWithControl:sLc];
        for (j = 0; j < nLoc; j++) {    // for 4 locs
            // (x, y) is original locacion
            // (x + xsft, y + ysft) is new location
            xsft[j] = *s;
            s += len;
            ysft[j] = *s;
            s += len;
            ctr = [self centerLocForPt:j view:[peIx current] nAng:nView];
            xi = ctr.x + xsft[j];
            yi = ctr.y + ysft[j];
            p[yi * xDim + xi] = mark;    // mark center of
        }
        [lc increment];
    }
// combine / filter and save
    pw_c = [self combinePWForLoop:ch withCoil:GE_Card_8_new];
    [pw_c pwFilter];
    [pw_c saveAsKOImage:@"pw_marker.img"];   // with marker, before correction
}

- (void)sideLobeFilter
{
    void    (^proc)(float *p, int xDim, int yDim);

    proc = ^(float *p, int xDim, int yDim) {
        int     i, j, ix;
        float   x, y;
        float   y0 = 0.66;
        float   w, wy, wy1, wy2;
        float   width = 0.15;

        for (i = ix = 0; i < yDim; i++) {
            y = ((float)i - yDim/2) * 2.0/yDim;
            wy1 = exp(-(y - y0) * (y - y0) / (width * width));
            wy2 = exp(-(y + y0) * (y + y0) / (width * width));
            wy = Rec_max(wy1, wy2);
            for (j = 0; j < xDim; j++, ix++) {
                x = ((float)j - xDim/2) * 2.0/xDim;
                w = wy * exp(-x * x / (width * width));
                p[ix] *= 1.0 - w;
            }
        }
    };
    [self apply2dProc:proc];
}

- (RecImage *)expandSliceBy:(int)ovr
{
    RecImage        *img;
    RecLoop         *yLp, *newY;     // PW: y is slice ...
    RecLoopControl  *srcLc, *dstLc;
    int             srcSkip, dstSkip;
    int             i, j, jj, k, n, len, len2;
    float           *p1, *p2;

    yLp = [self yLoop];
    len = [yLp dataLength];
    len2 = len + ovr*2;
    newY = [RecLoop loopWithDataLength:len2];

    img = [RecImage imageWithImage:self];
    [img replaceLoop:yLp withLoop:newY];
    srcSkip = [self skipSizeForLoop:yLp];
    dstSkip = [img skipSizeForLoop:newY];
    srcLc = [RecLoopControl controlForImage:self];
    dstLc = [RecLoopControl controlWithControl:srcLc forImage:img];
    [srcLc deactivateLoop:yLp];
    n = [srcLc loopLength];

    [srcLc rewind];
    [dstLc rewind];
    for (i = 0; i < n; i++) {
        p1 = [self currentDataWithControl:srcLc];
        p2 = [img currentDataWithControl:dstLc];
        for (k = 0; k < pixSize; k++) {
            for (j = 0; j < len; j++) {
                // center
                jj = j + ovr;
                p2[jj * dstSkip] = p1[j * srcSkip];
                // low
                jj = j + ovr - len;
                if (jj >= 0) {
                    p2[jj * dstSkip] = p1[j * srcSkip];
                }
                // hi
                jj = j + ovr + len;
                if (jj < len2) {
                    p2[jj * dstSkip] = p1[j * srcSkip];
                }
            }
            p1 += dataLength;
            p2 += [img dataLength];
        }
        [srcLc increment];
    }
    return img;
}

- (void)cropBy:(int)ovr
{
	RecLoop		*yLoop, *newYLoop;
	int			newYDim;

	yLoop = [self yLoop];
	newYDim = [yLoop dataLength] - ovr*2;
	newYLoop = [RecLoop loopWithDataLength:newYDim];

	[self replaceLoop:yLoop withLoop:newYLoop];
}

- (RecImage *)radCorrWith:(RecImage *)ref
{
    RecImage        *ref_r, *img_r;
    RecImage        *rFilt;
    RecLoop         *tLoop, *rLoop;
    RecLoopControl  *lc;
    NSMutableArray  *lpArray;
    int             dim;
    int             nTheta = 128;
    int             nRad = 8;
    int             i, j;
    float           w, *p;

// ### these images doesn't have common loops...
// make something like -(RecImage *)toPolarWithControl:(RecLoopControl *)lc;

//    img_r = [self toPolarWithNTheta:nTheta nRad:nRad];
//    ref_r = [ref toPolarWithNTheta:nTheta nRad:nRad];

// === make method for this part ====
    lpArray = [NSMutableArray arrayWithArray:[self loops]];
    dim = (int)[lpArray count];
    tLoop = [RecLoop loopWithDataLength:nTheta];
    rLoop = [RecLoop loopWithDataLength:nRad];
    [lpArray removeObjectAtIndex:dim-1];
    [lpArray removeObjectAtIndex:dim-2];
    [lpArray addObject:rLoop];
    [lpArray addObject:tLoop];

    img_r = [RecImage imageOfType:RECIMAGE_COMPLEX withLoopArray:lpArray];
    ref_r = [RecImage imageOfType:RECIMAGE_COMPLEX withLoopArray:lpArray];
// === make lib for thits part ====

    [img_r toPolarFrom:self];
    [ref_r toPolarFrom:ref];
    rFilt = [RecImage imageOfType:RECIMAGE_REAL withLoops:rLoop, tLoop, nil];
    lc = [RecLoopControl controlForImage:rFilt];
    [lc rewind];
    [lc deactivateInner];

// r filter    
    for (i = 0; i < nRad; i++) {
        w = (float)i/nRad;
        w = sin(M_PI * w);
        w *= w;
        p = [rFilt currentDataWithControl:lc];
        for (j = 0; j < nTheta; j++) {
            p[j] = w;
        }
        [lc increment];
    }
    [img_r scaleByImage:rFilt];
    [ref_r scaleByImage:rFilt];

//    return [img_r xyCorrelationWith:ref_r normalize:YES];
    return [img_r xyCorrelationWith:ref_r];
}

- (void)pwFilter
{
    BOOL    cpx = (type == RECIMAGE_COMPLEX);

    [self edgeFilter];
    if (!cpx) {
        [self makeComplex];
    }
    [self fft2d:REC_INVERSE];
    [self rhoFilter];
    [self sideLobeFilter];
    [self fft2d:REC_FORWARD];
    if (!cpx) {
        [self takeRealPart];
    }
}

- (void)edgeFilter
{
    void        (^proc)(float *p, int yDim, int skip);
    RecLoop     *yLoop = [self yLoop];
    float       *w, th;
    int         i, n = 8;

    w = (float *)malloc(sizeof(float) * n);
    for (i = 0; i < n; i++) {
        th = (float)i * 0.5 * M_PI / n;
        w[i] = sin(th);
    //    printf("%d %f\n", i, w[i]);
    }

    proc = ^(float *p, int yDim, int skip) {
        int     i, ix;
        float   m = (p[0] + p[(yDim - 1)* skip]) / 2.0;

        for (i = 0; i < n; i++) {
            ix = i * skip;
            p[ix] = (p[ix] - m) * w[i] + m;

            ix = (yDim - i - 1) * skip;
            p[ix] = (p[ix] - m) * w[i] + m;
        }
    };
    [self apply1dProc:proc forLoop:yLoop];
    free(w);
}

- (void)rhoFilter
{
    void    (^proc)(float *p, int xDim, int yDim);

    proc = ^(float *p, int xDim, int yDim) {
        int     i, j;
        float   rx, ry;

        for (i = 0; i < yDim; i++) {
            ry = (float)i - yDim/2;
            ry /= yDim;
            ry *= ry;
            for (j = 0; j < xDim; j++) {
                rx = (float)j - xDim/2;
                rx /= xDim;
                rx *= rx;
                p[i * xDim + j] *= sqrt(rx + ry);
            }
        }
    };
    [self apply2dProc:proc];
}

