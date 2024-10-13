# squid-approach
approach tuning + thermal
<h1>Nanonis API</h1>

Request message format:
Header:
40 bytes in size
command name (string) (32 bytes) name of the executed command 
Body size (int) (4 bytes) size of the message body in bytes
Send response back (unsigned int16) (2 bytes) defines if the server sends a message back to the client (=1) or not (=0)

Body:
contains argument values sent to setver

Response message format (iff response back 1):
Header40 bytes:
commandname (string)(32 bytes)
body size (int) (4 bytes)
Body: after the return arguments there is the error information:
Erorr status (unsigned int32)(4 bytes) returns 1=True if there is an error when executing the function
Error distribtion size (int)(4 bytes) returns the size of the error description
Error descritption (string): returns the description of the error
size depends on the function
## send,receive,Header functions

| Request message                                   ||||||
|---------------------------------------------------|----------------------------------------------------------------------------------|---------------------|------------|-----------------------|-----------| ||        Header                                 | Body                  |           |
| Command name                                      | Body size in bytes                                                               | Send response back  | Not Used   | Wait for newest data  |           |
| Size (Bytes)                                      | Fixed (32)                                                                       | Fixed (4)           | Fixed (2)  | Fixed (2)             | 4         |
| Readable value representation                     | FolMe.XYPosGet                                                                   | 4                   | True       | -                     | True      |
| Hex representation of string to be sent over TCP  | 466F 6C4D 652E 5859 506F 7347 6574 0000 0000 0000 0000 0000 0000 0000 0000 0000  | 0000 0004           | 0001       | 0000                  | 0000 0001 |


TODO:
check Header function, seems not exact
TODO:
change code to use tcpclient
seems more elegant