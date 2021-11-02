// mysql.cpp : Defines the entry point for the console application.
//
#include "string.h"
#include <iostream>
#include "C:\Program Files\MySQL\MySQL Server 8.0\include\mysql.h"
#pragma comment(lib,"libmySQL.lib")
using namespace std;
 MYSQL *conn;
 string mysqlsearchpromoter(string id)
 { 
	 string sss;
   conn=mysql_init(NULL); 
   if(!mysql_real_connect(conn, "localhost","root","123456","design",0,NULL,0))  
	{  
		cout<<"cannot connect datebase"; 
	    exit(1);
	 }  
//	 cout<<"连接成功";
	
	 char sql[100];
     sprintf(sql,"select sequence from promoter where id= \'%s\'",id.c_str());
	
	 MYSQL_RES *res_set; 
     MYSQL_ROW row; 
	  int r=mysql_query(conn,sql);
	  
	 
	  if(r==0)
	 {   
		 res_set = mysql_use_result(conn); 
		 row=mysql_fetch_row(res_set);
		 if(row==NULL){cout<<"ec"<<id.c_str()<<"is not in our database"; return ""; }
		 else{
			 sss=row[0];
		 }

	  }
	  mysql_close(conn);
	  return sss;

 }	
string mysqlsearchsd(string name)
 { 
	 string sss;
   conn=mysql_init(NULL); 
   if(!mysql_real_connect(conn, "localhost","root","123456","design",0,NULL,0))  
	{  
		cout<<"cannot connect datebase"; 
	    exit(1);
	 }  
//	 cout<<"连接成功";
	
	 char sql[100];
     sprintf(sql,"select seq from sd where name= \'%s\'",name.c_str());
	
	 MYSQL_RES *res_set; 
     MYSQL_ROW row; 
	  int r=mysql_query(conn,sql);
	  
	 
	  if(r==0)
	 {   
		 res_set = mysql_use_result(conn); 
		 row=mysql_fetch_row(res_set);
		 if(row==NULL){cout<<"ec"<<name.c_str()<<"is not in our database"; return ""; }
		 else{
			 sss=row[0];
		 }

	  }
	  mysql_close(conn);
	  return sss;

 }  

string mysqlsearchterminator(string name)
 { 
	 string sss;
   conn=mysql_init(NULL); 
   if(!mysql_real_connect(conn, "localhost","root","123456","design",3308,NULL,0))  {  
            cout<<"cannot connect datebase"; 
	    exit(1);
    }  
//	 cout<<"连接成功";
	
      char sql[100];
      sprintf(sql,"select sequence from terminater where name= \'%s\'",name.c_str());
	
	 MYSQL_RES *res_set; 
        MYSQL_ROW row; 
	  int r=mysql_query(conn,sql);
	  
	 
	  if(r==0)
	 {   
		 res_set = mysql_use_result(conn); 
		 row=mysql_fetch_row(res_set);
		 if(row==NULL){cout<<"ec"<<name.c_str()<<"is not in our database"; return ""; }
		 else{
			 sss=row[0];
		 }

	  }		
	  
	  mysql_close(conn);
	  return sss;

 }  

string mysqlres(string name)
 { 
	 string sss="";
   conn=mysql_init(NULL); 
   if(!mysql_real_connect(conn, "localhost","root","123456","design",3306,NULL,0))  {  
            cout<<"cannot connect datebase"; 
	    exit(1);
    }  
//	 cout<<"连接成功";
	
      char sql[100];
      sprintf(sql,"select seq from restriction where name= \'%s\'",name.c_str());
	
	 MYSQL_RES *res_set; 
        MYSQL_ROW row; 
	  int r=mysql_query(conn,sql);
	  
	 
	  if(r==0)
	 {   
		 res_set = mysql_use_result(conn); 
		 row=mysql_fetch_row(res_set);
		 if(row==NULL){
			 //使用second_name查一次
			 sprintf(sql, "select seq from restriction where second_name= \'%s\'", name.c_str());
			 r = mysql_query(conn, sql);
			 if (r == 0)
			 {
				 res_set = mysql_use_result(conn);
				 row = mysql_fetch_row(res_set);
				 if (row == NULL){ cout << "ec" << name.c_str() << "is not in our database"; return ""; }
				 else{
					 sss = row[0];
				 }

			 }
			 else{
				 cout << "ec" << name.c_str() << "is not in our database"; return "";
			 }
			 
		 }
	
		 else{
			 sss = row[0];
		 }

	  }
	  else{
		  cout << "ec" << name.c_str() << "is not in our database"; return "";
	  }
	  mysql_close(conn);
	  return sss;

 }  
