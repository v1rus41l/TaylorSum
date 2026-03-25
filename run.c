#include <stdio.h>
#include <math.h>

#define MAX_TERMS 1000000

typedef struct {
    float term;
    int n;
    float x;
} TaylorTerms;

void next_sin(TaylorTerms *s) {
    s->term *= -(s->x * s->x) / ((2.0f * s->n + 2) * (2.0f * s->n + 3));
    s->n++;
}

void next_cos(TaylorTerms *s) {
    s->term *= -(s->x * s->x) / ((2.0f * s->n + 1) * (2.0f * s->n + 2));
    s->n++;
}

void next_exp(TaylorTerms *s) {
    s->term *= s->x / (1.0f + s->n);
    s->n++;
}

void next_ln(TaylorTerms *s) {
    float t = (s->x - 1.0f) / (s->x + 1.0f);
    s->term *= (t * t) * (2.0f * s->n + 1) / (2.0f * s->n + 3);
    s->n++;
}

float forward_sum(TaylorTerms *s, void (*next)(TaylorTerms *)) {
    float sum = s->term;

    for (int i = 0; i < MAX_TERMS; i++) {
        next(s);
        sum += s->term;
    }

    return sum;
}

float back_sum (TaylorTerms *s, void (*next)(TaylorTerms *)) {
    float terms[MAX_TERMS];
    int index = 0;
    float sum = 0;

    terms[index] = s->term;
    index++;

    for (int i = 0; i < MAX_TERMS - 1; i++) {
        next(s);
        terms[index] = s->term;
        index++;
    }

    for (int i = index - 1; i >= 0; i--) {
         sum += terms[i];
    }

    return sum;
}

float kahan_algo(TaylorTerms *s, void (*next)(TaylorTerms *)) {
    float sum = s->term;
    float d = 0;

    for (int i = 0; i < MAX_TERMS; i++) {
        next(s);

        float a = s->term + d;
        float b = sum;

        float c = a + b;
        float t = c - b;

        d = a - t;
        sum = c;

    }
    return sum;
}

float period(float x, int *sign) {
    if (x < 0.0f) {
        *sign = -1;
        x = -x;
    } else {
        *sign = 1;
    }

    int count = (int)(x / (2.0f * (float)M_PI));
    x = x - count * 2.0f * (float)M_PI;

    return x;

}

int main() {
    float x;
    int sign;
    int func;
    int algo;
    float new_x;
    float sum;
    float from_math_value;

    printf("1 — sin(x)\n");
    printf("2 — cos(x)\n");
    printf("3 - exp(x)\n");
    printf("4 - ln(x)\n");
    printf("Choose function: ");
    scanf("%d", &func);

    printf("1 — direct\n");
    printf("2 — reverse\n");
    printf("3 — Kahan\n");
    printf("Choose method: ");
    scanf("%d", &algo);

    printf("Enter x: ");
    scanf("%f", &x);

    sign = 1;
    if (func == 1 || func == 2) {
        new_x = period(x, &sign);
    } else {
        new_x = x;
    }

    TaylorTerms s;
    s.n = 0;
    s.x = new_x;

    if (func == 1) {
        s.term = new_x;
    } else if (func == 4) {
        s.term = 2.0f * (x - 1.0f) / (x + 1.0f);
    } else {
        s.term = 1.0f;
    }

    void (*next)(TaylorTerms *);
    if (func == 1) {
        next = next_sin;
        from_math_value = sinf(x);
    } else if (func == 2) {
        next = next_cos;
        from_math_value = cosf(x);
    } else if (func == 3) {
        next = next_exp;
        from_math_value = expf(x);
    } else if (func == 4) {
        next = next_ln;
        from_math_value = logf(x);
    } else {
        printf("Error.\n");
        return 1;
    }

    if (algo == 1) {
        sum = forward_sum(&s, next);
    } else if (algo == 2) {
        sum = back_sum(&s, next);
    } else if (algo == 3) {
        sum = kahan_algo(&s, next);
    } else {
        printf("Error.\n");
        return 1;
    }

    if (func != 2) {
        sum *= sign;
    }

    printf("\n");
    printf("------------------------------\n");
    printf("Result: %.7f\n", sum);
    printf("Value from math.h: %.7f\n", from_math_value);
    printf("Error = %e\n", fabsf(sum - from_math_value));
    printf("------------------------------\n");
}