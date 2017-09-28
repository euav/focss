#ifndef COMMUNICATION_LINE_H_
#define COMMUNICATION_SYSTEM_H_

#include <ostream>
#include "focss/field.h"
#include "focss/component/fiber.h"
#include "focss/component/edfa.h"
#include "focss/solver/ssfm.h"

class CommunicationLine {
    unsigned long spans;
    SSFM ssfm;
    EDFA edfa;
    
    std::ostream* output_stream;

    
    
};

#endif // COMMUNICATION_LINE_H_
