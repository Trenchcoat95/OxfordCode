// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME EDepSimEventsProjectDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "EDepSimEventsProjectHeaders.h"

// Header files passed via #pragma extra_include

namespace {
  void TriggerDictionaryInitialization_EDepSimEventsProjectDict_Impl() {
    static const char* headers[] = {
"EDepSimEventsProjectHeaders.h",
0
    };
    static const char* includePaths[] = {
"/mnt/c/Linux/root6/root-bin/include",
"/mnt/c/Linux/root6/root-bin/etc",
"/mnt/c/Linux/root6/root-bin/etc/cling",
"/mnt/c/Linux/root6/root-bin/include",
"/mnt/c/Linux/GEANT/geant4_10_05_p01-install/include/Geant4",
"/mnt/c/Linux/Dune/edep-sim/io",
"/mnt/c/Linux/Dune/edep-sim/edep-gcc-7-x86_64-linux-gnu/io/",
"/mnt/c/Linux/root6/root-bin/include",
"/mnt/c/Users/fedeb/Documents/Universita/Federico 2020-2021/OxfordCode/EDepSimEvents/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "EDepSimEventsProjectDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "EDepSimEventsProjectDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "EDepSimEventsProjectHeaders.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("EDepSimEventsProjectDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_EDepSimEventsProjectDict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_EDepSimEventsProjectDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_EDepSimEventsProjectDict() {
  TriggerDictionaryInitialization_EDepSimEventsProjectDict_Impl();
}
