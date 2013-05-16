# $Id: queue.m4 $
# $Author: Magdalena Slawinska $
# $Date: 2010/06/18 $
#
# Autoconf macros defining the following substitution variables used by for NQS:
# subCMD - the name of queue (default is qsub)
# CLASS -  the name of class (default is qunlimitted)
#


AC_DEFUN([QUEUE_CMD],
[
	AC_ARG_WITH(queue,
	[  --with-queue          name of the queue used in NQS; default is qsub],
   	  user_queue=$withval,
	  user_queue="none")

if test ! x"$user_queue" = xnone; then
    subCMD="$user_queue" 
elif test ! x"$user_queue" = x ; then 
    subCMD="qsub"
else 
    subCMD="qsub"
fi

  AC_SUBST(subCMD)
  ])

AC_DEFUN([CLASS_CMD],
[
	AC_ARG_WITH(class,
	[  --with-class          name of the CLASS used in NQS default is qsub],
   	  user_CLASS=$withval,
	  user_CLASS="none")

if test ! x"$user_CLASS" = xnone; then
    CLASS="$user_CLASS" 
elif test ! x"$user_CLASS" = x ; then 
    CLASS="qunlimitted"
else 
    CLASS="qunlimitted"
fi

  AC_SUBST(CLASS)
  ])

#
# EOF
#
