//
//
//

// ===== old (not used) =======
- (RecImage *)expandSliceBy:(int)ovr;
- (void)cropBy:(int)ovr;
// filters ... not used anymore (along with 4pt method)
- (void)pwFilter;
- (void)edgeFilter;
- (void)rhoFilter;
- (void)sideLobeFilter;
// 4 pts correction scheme
- (RecImage *)detectShift4WithRef:(RecImage *)ref;
- (RecImage *)correctShift4:(RecImage *)sft mark:(BOOL)dbg over:(int)ovr;
- (void)markPoints:(RecImage *)sft;
- (NSPoint)centerLocForPt:(int)ix view:(int)pe nAng:(int)nView;
// radial correlation (for rotational correction)
- (RecImage *)radCorrWith:(RecImage *)ref;
- (void)toPolarFrom:(RecImage *)img;
