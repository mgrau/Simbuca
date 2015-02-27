#ifndef SMsgType_H
#define SMsgType_H

/// Message type enumeration
/**
 *          Written by Attila Krasznahorkay
 * Enumeration for classifying messages. Based on this classification,
 * SLogWriter can decide if/how the message should be shown.
 * (Naming taken from ATLAS offline.)
 */
enum SMsgType {
   VERBOSE = 1, /**< Type for the most detailed messages. Only for serious debugging. */
   DEBUG = 2,   /**< Type for debugging messages. A few messages per event allowed. */
   INFO = 3,    /**< Type for normal information messages. No messages in event processing! */
   WARNING = 4, /**< Type for smaller problems.  */
   ERROR = 5,   /**< Type for "real" problems. */
   FATAL = 6,   /**< Type for problems that should halt the execution. */
   ALWAYS = 7   /**< Type that should always be shown. */
};

#endif // SMsgType_H
