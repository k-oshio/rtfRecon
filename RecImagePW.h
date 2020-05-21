//
//  Paddle Wheel category
//

#import <RecKit/RecKit.h>

@interface	RecImage (PW)

// pre-processing (removing bad data)
- (void)remove2RR:(RecImage *)pw thres:(float)th;
- (void)removeDeepBreath:(RecImage *)sft thres:(float)th;
- (RecLoop *)sftCenter:(RecLoop *)lp;
- (RecLoop *)sftCenter:(RecLoop *)lp zeroFillTo:(int)zDim;
- (RecLoop *)halfFT2:(RecLoop *)lp;				// tmp... replace haflFT when done

// coil map
- (RecImage *)ftFitWithMask:(RecImage *)mask;

// gridding
- (void)radPhaseCorrWithGATab:(RecImage *)tab;	// in/out : raw [ch sl proj rd]
- (void)radPhaseCorr360;						// self-pcorr using 360 data (g.a. or not) 
- (RecImage *)reprojectWithMap:(RecImage *)radMap;    // img -> pw
- (RecImage *)pwToSin;

//- (RecImage *)goldenAngleTab;
- (void)initRadialTrajWithTab:(RecImage *)tab startPoint:(int)st;
- (void)initRadialMapWithTab:(RecImage *)tab;	// always full echo
- (RecImage *)mapForRadial;
// reformat
- (RecImage *)toCoronal;
- (RecImage *)toSagittal;

- (void)reorder:(RecLoop *)lp withTab:(RecImage *)tab;
- (void)invReorder:(RecLoop *)lp withTab:(RecImage *)tab;
- (void)initRadialTrajWithSortTab:(RecImage *)tab startPoint:(int)st;

// auto segmentation ---> fermi filter
- (void)pwWinWithRx:(float)rx ry:(float)ry d:(float)d x:(float)x y:(float)y;
- (void)imgWin;
- (void)SATx:(float)x w:(float)w th:(float)th;

// shift est in k-space for eddy current correction
// probably similar technique can be used for pw images
- (RecImage *)pDiffWith:(RecImage *)ref;
// shift correction using phase correlation
//- (RecImage *)corrToSftOld;	// current
- (RecImage *)corrToSft;	// top level (xy)
- (RecImage *)corrToSftZ;	// top level (z only)
- (RecImage *)corrToSftMC;	// multi-channel, xy		### not implemented yet
- (RecImage *)combineSftXY;	// svd or weighed sum		### not implemented yet
- (RecImage *)combineSftZ;	// svd ? imag only			### not implemented yet
- (RecImage *)detectShift;			// optical flow or laplacian
- (RecImage *)opticalFlow1d;
- (RecImage *)opticalFlow2dWithRef:(RecImage *)ref;
- (RecImage *)opticalFlow3dWithRef:(RecImage *)ref;
- (RecImage *)shiftFromK0;	// shift estimation from center-of-kspace view (POCS)
- (RecImage *)estShiftWithFrac:(float)fr slice:(int)slc nIter:(int)nIter err:(float *)err dbg:(BOOL)dbg;   // innter loop of above
- (RecImage *)stepCorrWithPW:(RecImage *)pw gridder:(RecGridder *)grid sft:(RecImage *)sft; // input(self) is img

// correction
- (RecImage *)correctShift:(RecImage *)sft;		// xy (Pw)
- (RecImage *)correctZShift:(RecImage *)sft;	// z

- (RecImage *)corr1dWithRef:(RecImage *)ref;
- (RecImage *)corrToSft1d:(RecLoop *)lp;
//- (RecImage *)correctShift1d:(RecImage *)sft forLoop:(RecLoop *)lp;


// stopping criterion
- (float)deltaWithPreviousSft:(RecImage *)sft;

// expand... increase size and fill with outermost value
- (void)expandLoop:(RecLoop *)lp to:(int)sz;

// testing...
//- (void)phaseUnwrap;		// discrete version
//- (RecImage *)unwrapEst;	// initial est using Laplacian
- (void)phaseUnwrap3D;
- (RecImage *)coilMap:(RecLoop *)ch;
- (void)echoFiltForLoop:(RecLoop *)lp;	// q&d filter for residual echo removal
- (void)ktPCA;		// testing...
- (void)svd_U:(RecImage **)uimg S:(RecImage **)simg Vt:(RecImage **)vtimg;

- (RecImage *)sortZSft;
- (RecImage *)binTabWithNBins:(int)nBins binSize:(int)sz;
- (RecImage *)breakLoop:(RecLoop *)lp intoBins:(RecImage *)bTab;
- (RecImage *)defVectorWithRef:(RecImage *)ref;
- (RecImage *)defToSftWithTab:(RecImage *)gaSort sft:(RecImage *)sft; // sft [blk, pe]
- (RecImage *)outerBlock:(RecBlock3)blk;
- (RecImage *)innerBlock:(RecBlock3)blk;
- (RecImage *)projectBlock:(RecBlock3)blk theta:(RecImage *)thTab;
- (RecImage *)selectSft:(RecImage *)sftIx;
- (RecImage *)selectSftF:(RecImage *)sftIx;
- (void)scaleVector:(RecVector)v shift:(RecImage *)sft rotate:(RecImage *)tab;

//- (RecImage *)loadEdge:(NSString *)path;

// warp
- (RecImage *)scaleMapToImage:(RecImage *)img;

// binned recon
- (void)scaleToBins;

// ==== testing... ======
// mag filter
- (void)magLPF:(float)w;
// polar filter
- (void)polarLPF:(float)w;

// combine with thres
- (RecImage *)combineForLoop:(RecLoop *)lp threshold:(float)th;

- (void)normalize2ForLoop:(RecLoop *)lp withMask:(RecImage *)msk;	// L2 norm
// ==== testing... ======

// debug
- (void)saveShift:(int)ix;
- (void)dumpShift;
- (void)saveComb:(RecLoop *)ch as:(NSString *)path;
- (void)printDefAt:(RecVector)v;
- (void)setMarker:(float)maxVal;

@end

// === c func
void		Rec_calcTransformParam(float *v, float *x, float *y, float *dx, float *dy, int nfree);
int			takeROI2d(RecImage *src, RecImage *dst, int xc, int yc, int slice);
RecImage *	makeSubImage2(RecImage *img, int xsz, int ysz, int xc, int yc);
RecImage *	makeSubImage3(RecImage *img, int xsz, int ysz, int zsz, int xc, int yc, int zc);
RecImage *	Rec_radDelayMap(RecImage *pw);
RecImage *	dispToMap2d(RecImage *disp);
RecImage *	dispToMap3d(RecImage *disp);
RecImage *	Rec_goldenAngleTab(RecLoop *pe, int nRef);
RecImage *	Rec_goldenAngleSort(RecImage *gaTab);
RecImage *	Rec_RadialTab(RecImage *tab);
RecImage *	Rec_load_edge(NSString *path);	// remove when below is done

float       Rec_find_nearest_peak(float *p, int skip, int len);

void		qsort_b_float(float *p, int n, int order);	// order 1:acsending, -1:descending
void        block_sort(RecImage *raw);	// not used anymore

