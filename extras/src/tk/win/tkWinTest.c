/*
 * tkWinTest.c --
 *
 *	Contains commands for platform specific tests for the Windows
 *	platform.
 *
 * Copyright (c) 1997 Sun Microsystems, Inc.
 * Copyright (c) 2000 by Scriptics Corporation.
 * Copyright (c) 2001 by ActiveState Corporation.
 *
 * See the file "license.terms" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#include "tkWinInt.h"

HWND tkWinCurrentDialog;

/*
 * Forward declarations of functions defined later in this file:
 */

static int		TestclipboardObjCmd(ClientData clientData,
			    Tcl_Interp *interp, int objc,
			    Tcl_Obj *const objv[]);
static int		TestwineventCmd(ClientData clientData,
			    Tcl_Interp *interp, int argc, const char **argv);
int			TkplatformtestInit(Tcl_Interp *interp);

/*
 *----------------------------------------------------------------------
 *
 * TkplatformtestInit --
 *
 *	Defines commands that test platform specific functionality for Windows
 *	platforms.
 *
 * Results:
 *	A standard Tcl result.
 *
 * Side effects:
 *	Defines new commands.
 *
 *----------------------------------------------------------------------
 */

int
TkplatformtestInit(
    Tcl_Interp *interp)		/* Interpreter to add commands to. */
{
    /*
     * Add commands for platform specific tests on MacOS here.
     */

    Tcl_CreateObjCommand(interp, "testclipboard", TestclipboardObjCmd,
	    (ClientData) Tk_MainWindow(interp), (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "testwinevent", TestwineventCmd,
	    (ClientData) Tk_MainWindow(interp), (Tcl_CmdDeleteProc *) NULL);
    return TCL_OK;
}

/*
 *----------------------------------------------------------------------
 *
 * AppendSystemError --
 *
 *	This routine formats a Windows system error message and places it into
 *	the interpreter result. Originally from tclWinReg.c.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

static void
AppendSystemError(
    Tcl_Interp *interp,		/* Current interpreter. */
    DWORD error)		/* Result code from error. */
{
    int length;
    WCHAR *wMsgPtr;
    const char *msg;
    char id[TCL_INTEGER_SPACE], msgBuf[24 + TCL_INTEGER_SPACE];
    Tcl_DString ds;
    Tcl_Obj *resultPtr = Tcl_GetObjResult(interp);

    if (Tcl_IsShared(resultPtr)) {
	resultPtr = Tcl_DuplicateObj(resultPtr);
    }
    length = FormatMessageW(FORMAT_MESSAGE_FROM_SYSTEM
	    | FORMAT_MESSAGE_ALLOCATE_BUFFER, NULL, error,
	    MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (WCHAR *) &wMsgPtr,
	    0, NULL);
    if (length == 0) {
	char *msgPtr;

	length = FormatMessageA(FORMAT_MESSAGE_FROM_SYSTEM
		| FORMAT_MESSAGE_ALLOCATE_BUFFER, NULL, error,
		MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (char *) &msgPtr,
		0, NULL);
	if (length > 0) {
	    wMsgPtr = (WCHAR *) LocalAlloc(LPTR, (length + 1) * sizeof(WCHAR));
	    MultiByteToWideChar(CP_ACP, 0, msgPtr, length + 1, wMsgPtr,
		    length + 1);
	    LocalFree(msgPtr);
	}
    }
    if (length == 0) {
	if (error == ERROR_CALL_NOT_IMPLEMENTED) {
	    strcpy(msgBuf, "function not supported under Win32s");
	} else {
	    sprintf(msgBuf, "unknown error: %ld", error);
	}
	msg = msgBuf;
    } else {
	Tcl_Encoding encoding;
	char *msgPtr;

	encoding = Tcl_GetEncoding(NULL, "unicode");
	Tcl_ExternalToUtfDString(encoding, (char *) wMsgPtr, -1, &ds);
	Tcl_FreeEncoding(encoding);
	LocalFree(wMsgPtr);

	msgPtr = Tcl_DStringValue(&ds);
	length = Tcl_DStringLength(&ds);

	/*
	 * Trim the trailing CR/LF from the system message.
	 */

	if (msgPtr[length-1] == '\n') {
	    --length;
	}
	if (msgPtr[length-1] == '\r') {
	    --length;
	}
	msgPtr[length] = 0;
	msg = msgPtr;
    }

    sprintf(id, "%ld", error);
    Tcl_SetErrorCode(interp, "WINDOWS", id, msg, (char *) NULL);
    Tcl_AppendToObj(resultPtr, msg, length);
    Tcl_SetObjResult(interp, resultPtr);

    if (length != 0) {
	Tcl_DStringFree(&ds);
    }
}

/*
 *----------------------------------------------------------------------
 *
 * TestclipboardObjCmd --
 *
 *	This function implements the testclipboard command. It provides a way
 *	to determine the actual contents of the Windows clipboard.
 *
 * Results:
 *	A standard Tcl result.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

static int
TestclipboardObjCmd(
    ClientData clientData,	/* Main window for application. */
    Tcl_Interp *interp,		/* Current interpreter. */
    int objc,			/* Number of arguments. */
    Tcl_Obj *const objv[])	/* Argument values. */
{
    HGLOBAL handle;
    char *data;
    int code = TCL_OK;

    if (objc != 1) {
	Tcl_WrongNumArgs(interp, 1, objv, (char *) NULL);
	return TCL_ERROR;
    }
    if (OpenClipboard(NULL)) {
	/*
	 * We could consider using CF_UNICODETEXT on NT, but then we
	 * would have to convert it from External.  Instead we'll just
	 * take this and do "bytestring" at the Tcl level for Unicode
	 * inclusive text
	 */
	handle = GetClipboardData(CF_TEXT);
	if (handle != NULL) {
	    data = GlobalLock(handle);
	    Tcl_AppendResult(interp, data, (char *) NULL);
	    GlobalUnlock(handle);
	} else {
	    Tcl_AppendResult(interp, "null clipboard handle", (char *) NULL);
	    code = TCL_ERROR;
	}
	CloseClipboard();
	return code;
    } else {
	Tcl_AppendResult(interp, "couldn't open clipboard: ", (char *) NULL);
	AppendSystemError(interp, GetLastError());
	return TCL_ERROR;
    }
    return TCL_OK;
}

/*
 *----------------------------------------------------------------------
 *
 * TestwineventCmd --
 *
 *	This function implements the testwinevent command. It provides a way
 *	to send messages to windows dialogs.
 *
 * Results:
 *	A standard Tcl result.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

static int
TestwineventCmd(
    ClientData clientData,	/* Main window for application. */
    Tcl_Interp *interp,		/* Current interpreter. */
    int argc,			/* Number of arguments. */
    const char **argv)		/* Argument strings. */
{
    HWND hwnd = 0;
    HWND child = 0;
    int id;
    char *rest;
    UINT message;
    WPARAM wParam;
    LPARAM lParam;
    static const TkStateMap messageMap[] = {
	{WM_LBUTTONDOWN,	"WM_LBUTTONDOWN"},
	{WM_LBUTTONUP,		"WM_LBUTTONUP"},
	{WM_CHAR,		"WM_CHAR"},
	{WM_GETTEXT,		"WM_GETTEXT"},
	{WM_SETTEXT,		"WM_SETTEXT"},
	{WM_COMMAND,            "WM_COMMAND"},
	{-1,			NULL}
    };

    if ((argc == 3) && (strcmp(argv[1], "debug") == 0)) {
	int b;

	if (Tcl_GetBoolean(interp, argv[2], &b) != TCL_OK) {
	    return TCL_ERROR;
	}
	TkWinDialogDebug(b);
	return TCL_OK;
    }

    if (argc < 4) {
	return TCL_ERROR;
    }

    hwnd = (HWND) strtol(argv[1], &rest, 0);
    if (rest == argv[1]) {
	hwnd = FindWindow(NULL, argv[1]);
	if (hwnd == NULL) {
	    Tcl_SetResult(interp, "no such window", TCL_STATIC);
	    return TCL_ERROR;
	}
    }
    UpdateWindow(hwnd);

    id = strtol(argv[2], &rest, 0);
    if (rest == argv[2]) {
	char buf[256];

	child = GetWindow(hwnd, GW_CHILD);
	while (child != NULL) {
	    SendMessage(child, WM_GETTEXT, (WPARAM) sizeof(buf), (LPARAM) buf);
	    if (strcasecmp(buf, argv[2]) == 0) {
		id = GetDlgCtrlID(child);
		break;
	    }
	    child = GetWindow(child, GW_HWNDNEXT);
	}
	if (child == NULL) {
	    Tcl_AppendResult(interp, "could not find a control matching \"",
		argv[2], "\"", NULL);
	    return TCL_ERROR;
	}
    }
    message = TkFindStateNum(NULL, NULL, messageMap, argv[3]);
    if (message < 0) {
	message = strtol(argv[3], NULL, 0);
    }
    wParam = 0;
    lParam = 0;

    if (argc > 4) {
	wParam = strtol(argv[4], NULL, 0);
    }
    if (argc > 5) {
	lParam = strtol(argv[5], NULL, 0);
    }

    switch (message) {
	case WM_GETTEXT: {
	    Tcl_DString ds;
	    char buf[256];

	    GetDlgItemText(hwnd, id, buf, 256);
	    Tcl_ExternalToUtfDString(NULL, buf, -1, &ds);
	    Tcl_AppendResult(interp, Tcl_DStringValue(&ds), NULL);
	    Tcl_DStringFree(&ds);
	    break;
	}
	case WM_SETTEXT: {
	    Tcl_DString ds;
	    BOOL result;

	    Tcl_UtfToExternalDString(NULL, argv[4], -1, &ds);
	    result = SetDlgItemText(hwnd, id, Tcl_DStringValue(&ds));
	    Tcl_DStringFree(&ds);
	    if (result == 0) {
		Tcl_SetResult(interp, "failed to send text to dialog: ", TCL_STATIC);
		AppendSystemError(interp, GetLastError());
		return TCL_ERROR;
	    }
	    break;
	}
	case WM_COMMAND: {
	    char buf[TCL_INTEGER_SPACE];
	    if (argc < 5) {
		wParam = MAKEWPARAM(id, 0);
		lParam = (LPARAM)child;
	    }
	    sprintf(buf, "%d", (int) SendMessageA(hwnd, message, wParam, lParam));
	    Tcl_SetResult(interp, buf, TCL_VOLATILE);
	    break;
	}
	default: {
	    char buf[TCL_INTEGER_SPACE];

	    sprintf(buf, "%d",
		    (int) SendDlgItemMessage(hwnd, id, message, wParam, lParam));
	    Tcl_SetResult(interp, buf, TCL_VOLATILE);
	    break;
	}
    }
    return TCL_OK;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 4
 * fill-column: 78
 * End:
 */
