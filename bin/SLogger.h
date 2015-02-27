#ifndef SLogger_H
#define SLogger_H

// STL include(s):
#include <string>
#include <sstream>

// Local include(s):
#include "SMsgType.h"
#include "SLogWriter.h"

#ifdef __GUI_ON__
  #include <QObject>
  #include <QMessageBox>
  #include <QApplication>
#endif


// Forward declaration(s):
class SLogWriter;

/**
 *   @short Universal message logging class
 *
 *          Written by Attila Krasznahorkay
 *          Class that should be used in the whole framework for printing
 *          messages on the terminal. It produces nicely formatted log
 *          messages using the SLogWriter class.
 *
 *          It prints the source of each message which makes reading
 *          debugging messages a bit easier. Every message has a type.
 *          This type controls how/if they should be printed to the
 *          terminal with the current configuration.
 *
 *     @see SLogWriter
 *     @see SMsgType
 * @version $Revision: 344 $
 */
class SLogger : public std::ostringstream {

public:
   /// Constructor with a name of the parent object
   SLogger( const std::string& source );
   /// Copy constructor
   SLogger( const SLogger& parent );

   /// Set the source name using a simple string
   void SetSource( const std::string& source );

   /// Get the source string of the logger
   const char* GetSource() const;

   /// Copy operator
   SLogger& operator= ( const SLogger& parent );

   /// Stream modifier to send a message
   static SLogger& endmsg( SLogger& logger );

   /// Operator accepting SLogger stream modifiers
   SLogger& operator<< ( SLogger& ( *_f )( SLogger& ) );
   /// Operator accepting std::ostream stream modifiers
   SLogger& operator<< ( std::ostream& ( *_f )( std::ostream& ) );
   /// Operator accepting std::ios stream modifiers
   SLogger& operator<< ( std::ios& ( *_f )( std::ios& ) );

   /// Operator accepting message type setting
   SLogger& operator<< ( SMsgType type );

   /// Operator accepting basically any kind of argument
   /**
    * SLogger was designed to give all the features that std::ostream
    * objects usually provide. This operator handles all kinds of
    * arguments and passes it on to the std::ostringstream base class.
    *
    * @param arg Any type of argument that std::ostream can handle
    * @returns This same object
    */
   template < class T > SLogger& operator<< ( T arg ) {
      if( m_activeType >= m_logWriter->GetMinType() ) {
         ( * ( std::ostringstream* ) this ) << arg;
      }
      return *this;
   }

   /// Old style message sender function
   void Send( SMsgType type, const std::string& message ) const;

private:
   /// Internal function for sending the message to the console
   void Send();

   std::string    m_strSource; ///< Name of the object sending the messages
   SLogWriter*    m_logWriter; ///< Pointer to the log writer object
   SMsgType       m_activeType; ///< Currently active message type

}; // class SLogger

//////////////////////////////////////////////////////////////////////
//                                                                  //
//   To speed up the code a bit, the following operators are        //
//   declared 'inline'.                                             //
//                                                                  //
//////////////////////////////////////////////////////////////////////

/**
 * This operator handles all stream modifiers that have been written
 * to work on SLogger objects specifically. Right now there is basically
 * only the SLogger::endmsg stream modifier that is such.
 */
inline SLogger& SLogger::operator<< ( SLogger& ( *_f )( SLogger& ) ) {

   return ( _f )( *this );
}

/**
 * This operator handles all stream modifiers that have been written
 * to work on std::ostream objects. Most of the message formatting
 * modifiers are such.
 */
inline SLogger& SLogger::operator<< ( std::ostream& ( *_f )( std::ostream& ) ) {

   if( m_activeType >= m_logWriter->GetMinType() ) {
      ( _f )( *this );
   }
   return *this;
}

/**
 * This operator handles all stream modifiers that have been written
 * to work on std::ios objects. I have to admit I don't remember exactly
 * which operators these are, but some formatting operations need this.
 */
inline SLogger& SLogger::operator<< ( std::ios& ( *_f )( std::ios& ) ) {

   if( m_activeType >= m_logWriter->GetMinType() ) {
      ( _f )( *this );
   }
   return *this;
}

/**
 * Messages have a type, defined by the SMsgType enumeration. This operator
 * allows the user to write intuitive message lines in the code like this:
 *
 * <code>
 *   logger << INFO << "This is an info message" << SLogger::endmsg;
 * </code>
 */
inline SLogger& SLogger::operator<< ( SMsgType type ) {

   m_activeType = type;
   return *this;
}

// This is a GCC extension for getting the name of the current function.
#if defined( __GNUC__ )
#   define SLOGGER_FNAME __PRETTY_FUNCTION__
#else
#   define SLOGGER_FNAME ""
#endif

/// Common prefix for the non-usual messages
/**
 * The idea is that a regular user usually only wants to see DEBUG, INFO
 * and some WARNING messages. So those should be reasonably short. On the other
 * hand serious warnings (ERROR, FATAL) or VERBOSE messages should be as precise
 * as possible.
 *
 * So I stole the idea from Athena (what a surprise...) to have a few macros
 * which produce messages with a common formatting. This macro provides the
 * prefix for all the messages.
 */
#define SLOGGER_REPORT_PREFIX \
   __FILE__ << ":" << __LINE__ << " (" << SLOGGER_FNAME << "): "

/// Convenience macro for reporting VERBOSE messages in the code
/**
 * This macro is very similar to the REPORT_MESSAGE macros of Athena. It prints
 * a nicely formatted output that specifies both the exact function name where
 * the message was printed, and also the filename:line combination. It can be
 * used like a regular function inside cycles:
 *
 * <code>
 *   REPORT_VERBOSE( "This is a verbose message with a number: " << number );
 * </code>
 */
#define REPORT_VERBOSE( MESSAGE ) \
   m_logger << VERBOSE << SLOGGER_REPORT_PREFIX << MESSAGE << SLogger::endmsg

/// Convenience macro for reporting ERROR messages in the code
/**
 * This macro is very similar to the REPORT_MESSAGE macros of Athena. It prints
 * a nicely formatted output that specifies both the exact function name where
 * the message was printed, and also the filename:line combination. It can be
 * used like a regular function inside cycles:
 *
 * <code>
 *   REPORT_ERROR( "A serious error message" );
 * </code>
 */
#define REPORT_ERROR( MESSAGE ) \
   m_logger << ERROR << SLOGGER_REPORT_PREFIX << MESSAGE << SLogger::endmsg

/// Convenience macro for reporting FATAL messages in the code
/**
 * This macro is very similar to the REPORT_MESSAGE macros of Athena. It prints
 * a nicely formatted output that specifies both the exact function name where
 * the message was printed, and also the filename:line combination. It can be
 * used like a regular function inside cycles:
 *
 * <code>
 *   REPORT_FATAL( "A very serious error message" );
 * </code>
 */
#define REPORT_FATAL( MESSAGE ) \
   m_logger << FATAL << SLOGGER_REPORT_PREFIX << MESSAGE << SLogger::endmsg

#endif // SLogger_H
