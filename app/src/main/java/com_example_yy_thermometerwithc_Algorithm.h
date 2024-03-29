/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
/* Header for class com_example_yy_thermometerwithc_Algorithm */

#ifndef _Included_com_example_yy_thermometerwithc_Algorithm
#define _Included_com_example_yy_thermometerwithc_Algorithm
#ifdef __cplusplus
extern "C" {
#endif
/*
 * Class:     com_example_yy_thermometerwithc_Algorithm
 * Method:    correlationJni
 * Signature: (I[D[D)I
 */
JNIEXPORT jint JNICALL Java_com_example_yy_thermometerwithc_Algorithm_correlationJni
  (JNIEnv *, jclass, jint, jdoubleArray, jdoubleArray);

/*
 * Class:     com_example_yy_thermometerwithc_Algorithm
 * Method:    mixFrequenceJni
 * Signature: ([D[D)[D
 */
JNIEXPORT jdoubleArray JNICALL Java_com_example_yy_thermometerwithc_Algorithm_mixFrequenceJni
  (JNIEnv *, jclass, jdoubleArray, jdoubleArray);

/*
 * Class:     com_example_yy_thermometerwithc_Algorithm
 * Method:    fftJni
 * Signature: ([D)[D
 */
JNIEXPORT jdoubleArray JNICALL Java_com_example_yy_thermometerwithc_Algorithm_fftJni
  (JNIEnv *, jclass, jdoubleArray);

/*
 * Class:     com_example_yy_thermometerwithc_Algorithm
 * Method:    maxIndexInRangeJni
 * Signature: ([DII)I
 */
JNIEXPORT jint JNICALL Java_com_example_yy_thermometerwithc_Algorithm_maxIndexInRangeJni
  (JNIEnv *, jclass, jdoubleArray, jint, jint);

/*
 * Class:     com_example_yy_thermometerwithc_Algorithm
 * Method:    firHelperJni
 * Signature: ([D[D[S)[D
 */
JNIEXPORT jdoubleArray JNICALL Java_com_example_yy_thermometerwithc_Algorithm_firHelperJni
  (JNIEnv *, jclass, jdoubleArray, jdoubleArray, jshortArray);

/*
 * Class:     com_example_yy_thermometerwithc_Algorithm
 * Method:    normolizeArray
 * Signature: ([S)[D
 */
JNIEXPORT jdoubleArray JNICALL Java_com_example_yy_thermometerwithc_Algorithm_normolizeArray
  (JNIEnv *, jclass, jshortArray);

#ifdef __cplusplus
}
#endif
#endif
