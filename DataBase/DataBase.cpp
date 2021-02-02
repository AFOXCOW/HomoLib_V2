// mysql.cpp : Defines the entry point for the console application.
#include "DataBase.h"
MYSQL* conn;

string mysqlres(string ResSite_Name)
{
	string ResSite = "";
	conn = mysql_init(NULL);
	if (!mysql_real_connect(conn, "localhost", "root", "123456", "design", 3306, NULL, 0)) {
		cout << "cannot connect datebase"<<endl;
		exit(1);
	}
	char sql[100];
	sprintf_s(sql,sizeof(sql),"select seq from restriction where name= \'%s\'", ResSite_Name.c_str());
	MYSQL_RES* res_set;
	MYSQL_ROW row;
	int Failed = mysql_query(conn, sql);
	if (!Failed){
		res_set = mysql_use_result(conn);
		row = mysql_fetch_row(res_set);
		if (row == NULL) {
			cout << "Restriction Name: " << ResSite_Name.c_str() << " is not in our database";
			return "";
		}
		else {
			ResSite = row[0];
		}
	}
	else {
		cout << "Restriction Name: " << ResSite_Name.c_str() << " is not in our database"; 
		return "";
	}
	mysql_close(conn);
	return ResSite;
}




























/* Not Used

string mysqlres(string name)
{
	string sss = "";
	conn = mysql_init(NULL);
	if (!mysql_real_connect(conn, "localhost", "root", "123456", "design", 3306, NULL, 0)) {
		cout << "cannot connect datebase";
		exit(1);
	}
	char sql[100];
	sprintf(sql, "select seq from restriction where name= \'%s\'", name.c_str());
	MYSQL_RES* res_set;
	MYSQL_ROW row;
	int r = mysql_query(conn, sql);
	if (r == 0)
	{
		res_set = mysql_use_result(conn);
		row = mysql_fetch_row(res_set);
		if (row == NULL) {
			//使用second_name查一次
			sprintf(sql, "select seq from restriction where second_name= \'%s\'", name.c_str());
			r = mysql_query(conn, sql);
			if (r == 0)
			{
				res_set = mysql_use_result(conn);
				row = mysql_fetch_row(res_set);
				if (row == NULL) { cout << "ec" << name.c_str() << "is not in our database"; return ""; }
				else {
					sss = row[0];
				}
			}
			else {
				cout << "ec" << name.c_str() << "is not in our database"; return "";
			}
		}
		else {
			sss = row[0];
		}
	}
	else {
		cout << "ec" << name.c_str() << "is not in our database"; return "";
	}
	mysql_close(conn);
	return sss;
}


string mysqlsearchpromoter(string id)
{
	string sss;
	conn = mysql_init(NULL);
	if (!mysql_real_connect(conn, "localhost", "root", "123456", "design", 0, NULL, 0))
	{
		cout << "cannot connect datebase";
		exit(1);
	}
	char sql[100];
	sprintf(sql, "select sequence from promoter where id= \'%s\'", id.c_str());
	MYSQL_RES* res_set;
	MYSQL_ROW row;
	int r = mysql_query(conn, sql);
	if (r == 0)
	{
		res_set = mysql_use_result(conn);
		row = mysql_fetch_row(res_set);
		if (row == NULL) {
			cout << "ec" << id.c_str() << "is not in our database"; 
			return ""; 
		}
		else {
			sss = row[0];
		}
	}
	mysql_close(conn);
	return sss;

}
string mysqlsearchsd(string name)
{
	string sss;
	conn = mysql_init(NULL);
	if (!mysql_real_connect(conn, "localhost", "root", "123456", "design", 0, NULL, 0))
	{
		cout << "cannot connect datebase";
		exit(1);
	}
	char sql[100];
	sprintf(sql, "select seq from sd where name= \'%s\'", name.c_str());
	MYSQL_RES* res_set;
	MYSQL_ROW row;
	int r = mysql_query(conn, sql);
	if (r == 0)
	{
		res_set = mysql_use_result(conn);
		row = mysql_fetch_row(res_set);
		if (row == NULL) { cout << "ec" << name.c_str() << "is not in our database"; return ""; }
		else {
			sss = row[0];
		}
	}
	mysql_close(conn);
	return sss;
}

string mysqlsearchterminator(string name)
{
	string sss;
	conn = mysql_init(NULL);
	if (!mysql_real_connect(conn, "localhost", "root", "123456", "design", 3308, NULL, 0)) {
		cout << "cannot connect datebase";
		exit(1);
	}
	char sql[100];
	sprintf(sql, "select sequence from terminater where name= \'%s\'", name.c_str());
	MYSQL_RES* res_set;
	MYSQL_ROW row;
	int r = mysql_query(conn, sql);
	if (r == 0)
	{
		res_set = mysql_use_result(conn);
		row = mysql_fetch_row(res_set);
		if (row == NULL) { cout << "ec" << name.c_str() << "is not in our database"; return ""; }
		else {
			sss = row[0];
		}
	}
	mysql_close(conn);
	return sss;
}

*/