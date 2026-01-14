//
// Calculate biquad coefficients
//
struct biquad_coeff {
	float b0, b1, b2;
	float a1, a2;
};

struct biquad_state {
	float x[2], y[2];
};

struct biquad {
	struct biquad_coeff coeff;
	struct biquad_state state;
};

// Direct form 1 may need more state than the "canonical" DF2,
// but gets noisy much quicker.
static inline float biquad_step_df1(struct biquad_coeff *c, float in, float x[2], float y[2])
{
	float out = c->b0*in + c->b1*x[0] + c->b2*x[1] - c->a1*y[0] - c->a2*y[1];
	x[1] = x[0]; x[0] = in;
	y[1] = y[0]; y[0] = out;
	return out;
}

static inline float _biquad_step(struct biquad_coeff *c, struct biquad_state *s, float x0)
{
	return biquad_step_df1(c, x0, s->x, s->y);
}

static inline void _biquad_lpf(struct biquad_coeff *res, float f, float Q)
{
	struct sincos w0 = fastsincos(f/SAMPLES_PER_SEC);
	float alpha = w0.sin/(2*Q);
	float a0_inv = 1/(1 + alpha);
	float b1 = (1 - w0.cos) * a0_inv;

	res->b0 = b1 / 2;
	res->b1 = b1;
	res->b2 = b1 / 2;
	res->a1 = -2*w0.cos	* a0_inv;
	res->a2 = (1 - alpha)	* a0_inv;
}

static inline void _biquad_hpf(struct biquad_coeff *res, float f, float Q)
{
	struct sincos w0 = fastsincos(f/SAMPLES_PER_SEC);
	float alpha = w0.sin/(2*Q);
	float a0_inv = 1/(1 + alpha);
	float b1 = (1 + w0.cos) * a0_inv;

	res->b0 = b1 / 2;
	res->b1 = -b1;
	res->b2 = b1 / 2;
	res->a1 = -2*w0.cos	* a0_inv;
	res->a2 = (1 - alpha)	* a0_inv;
}

static inline void _biquad_notch_filter(struct biquad_coeff *res, float f, float Q)
{
	struct sincos w0 = fastsincos(f/SAMPLES_PER_SEC);
	float alpha = w0.sin/(2*Q);
	float a0_inv = 1/(1 + alpha);

	res->b0 = 1 		* a0_inv;
	res->b1 = -2*w0.cos	* a0_inv;
	res->b2 = 1		* a0_inv;
	res->a1 = -2*w0.cos	* a0_inv;
	res->a2 = (1 - alpha)	* a0_inv;
}

static inline void _biquad_bpf_peak(struct biquad_coeff *res, float f, float Q)
{
	struct sincos w0 = fastsincos(f/SAMPLES_PER_SEC);
	float alpha = w0.sin/(2*Q);
	float a0_inv = 1/(1 + alpha);

	res->b0 = Q*alpha	* a0_inv;
	res->b1 = 0;
	res->b2 = -Q*alpha	* a0_inv;
	res->a1 = -2*w0.cos	* a0_inv;
	res->a2 = (1 - alpha)	* a0_inv;
}

static inline void _biquad_bpf(struct biquad_coeff *res, float f, float Q)
{
	struct sincos w0 = fastsincos(f/SAMPLES_PER_SEC);
	float alpha = w0.sin/(2*Q);
	float a0_inv = 1/(1 + alpha);

	res->b0 = alpha		* a0_inv;
	res->b1 = 0;
	res->b2 = -alpha	* a0_inv;
	res->a1 = -2*w0.cos	* a0_inv;
	res->a2 = (1 - alpha)	* a0_inv;
}

static inline void _biquad_allpass_filter(struct biquad_coeff *res, float f, float Q)
{
	struct sincos w0 = fastsincos(f/SAMPLES_PER_SEC);
	float alpha = w0.sin/(2*Q);
	float a0_inv = 1/(1 + alpha);

	res->b0 = (1 - alpha)	* a0_inv;
	res->b1 = (-2*w0.cos)	* a0_inv;
	res->b2 = 1;		// Same as a0
	res->a1 = res->b1;
	res->a2 = res->b0;
}

static inline float biquad_step(struct biquad *bq, float x0)
{ return _biquad_step(&bq->coeff, &bq->state, x0); }

#define biquad_lpf(bq,f,Q) _biquad_lpf(&(bq)->coeff,f,Q)
#define biquad_hpf(bq,f,Q) _biquad_hpf(&(bq)->coeff,f,Q)
#define biquad_notch_filter(bq,f,Q) _biquad_notch_filter(&(bq)->coeff,f,Q)
#define biquad_bpf_peak(bq,f,Q) _biquad_bpf_peak(&(bq)->coeff,f,Q)
#define biquad_bpf(bq,f,Q) _biquad_bpf(&(bq)->coeff,f,Q)
#define biquad_allpass_filter(bq,f,Q) _biquad_allpass_filter(&(bq)->coeff,f,Q)
