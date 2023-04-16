using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Reflection.Emit;
using System.Text;
using System.Windows.Forms;
using static System.Windows.Forms.VisualStyles.VisualStyleElement.Button;


namespace LAB_CHM_2023_3_1
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private void тестоваяЗадачаToolStripMenuItem_Click(object sender, EventArgs e)
        {

        }

        // =========  Functions =========
        double u1(double x, double y) // U* Решение тестовой задачи
        {
            return Math.Exp(1 - Math.Pow(x, 2) - Math.Pow(y, 2));
        }
        double f1(double x, double y) // Функция полученная через Лапласса
        {
            return (2*(2*Math.Pow(x,2)-1)*Math.Exp(0-Math.Pow(x,2)-Math.Pow(y,2)+1)+ 2 * (2 * Math.Pow(y, 2) - 1) * Math.Exp(0 - Math.Pow(x, 2) - Math.Pow(y, 2) + 1));
        }

        double f2(double x, double y) // F*
        {
            return Math.Abs(Math.Pow(x,2) - Math.Pow(y, 2));
        }

        double mu1(double y) //Граничное условие 1
        {
            return -y*y+1;
        }

        double mu2(double y) //Граничное условие 2
        {
            return (1-y*y)*Math.Exp(y);
        }

        double mu3(double x) //Граничное условие 3
        {
            return 1-x*x;
        }

        double mu4(double x) //Граничное условие 4
        {
            return 1 - x * x;
        }

        private void основнаяЗадачаToolStripMenuItem_Click(object sender, EventArgs e)
        {

        }

        private void label46_Click(object sender, EventArgs e)
        {

        }

        private void groupBox2_Enter(object sender, EventArgs e)
        {

        }

        private void button1_Click(object sender, EventArgs e)
        {
            int n = Convert.ToInt32(textBox1.Text);
            int m = Convert.ToInt32(textBox2.Text);
            int N_max = Convert.ToInt32(textBox3.Text);
            double Eps = Convert.ToDouble(textBox4.Text);
            double h = 2.0 / (double)n, k = 2.0 / (double)m; //Шаги по x и y
            double h2 = 1.0 / (h * h), k2 = 1.0 / (k * k);
            double A = -2 * (h2 + k2);
            double[][] v;
            double[][] f;
            double[][] u;
            double[][] hv;
            double[][] R;
            double[] x, y;
            int p = 0; //Текущее число итераций
            char[] buffer = new char [100];
            double MaxPogr = 0.0;
            double Pogr;
            double MaxF = 0.0;
            double maxR1 = 0.0;

            x = new double[n + 1];
            y = new double[m + 1];
            v = new double[n + 1][];
            hv = new double[n + 1][];
            R = new double[n + 1][];
            f = new double[n + 1][];
            u = new double[n + 1][];

            for (int i = 0; i <= n; i++)
            {
                v[i] = new double[m + 1];
                f[i] = new double[m + 1];
                u[i] = new double[m + 1];
                hv[i] = new double[m + 1];
                R[i] = new double[m + 1];
            }
            
            for (int i = 0; i <= n; i++)  //Заполнение массива x
            {
                x[i] = -1 + i * h;
              
            }
     
            for (int j = 0; j <= m; j++)  //Заполнение массива y
            {
                y[j] = -1 + j * k;// БЫЛО 2+j
            
            }



            for (int j = 0; j <= m; j++)            //Заполнение массивов f и u
            {
                for (int i = 0; i <= n; i++)
                {
                    f[i][j] = f1(x[i], y[j]);
                    u[i][j] = u1(x[i], y[j]);
                    hv[i][j] = -f[i][j];
                    MaxF += f[i][j] * f[i][j];
                }
            }

            MaxF = Math.Sqrt(MaxF);

            for (int j = 0; j <= m; j++)  //Заполнение граничных условий в массив v
            {
                v[0][j] = u1(-1, y[j]);
                v[n][j] = u1(1, y[j]);
                R[0][j] = 0.0;
                R[n][j] = 0.0;
            }

            for (int i = 0; i <= n; i++)  //Заполнение граничных условий в массив v
            {
                v[i][0] = u1(x[i], -1);
                v[i][m] = u1(x[i], 1);
                R[i][0] = 0.0;
                R[i][m] = 0.0;
            }

            for (int j = 1; j < m; j++)    //Нулевое начальное приближение
            {
                for (int i = 1; i < n; i++)
                {
                    v[i][j] = 0.0;
                }
            }

            // UpRelaxMethod
            double temp, prev, currentEps;
            double Eps_max;
            double w = 1.99;
            while (true)
            {
                Eps_max = 0.0;
                for (int j = 1; j < m; j++)
                {
                    for (int i = 1; i < n; i++)
                    {
                        prev = v[j][i];
                        temp = A * prev + h2 * (v[j][i - 1] + v[j][i + 1]) + k2 * (v[j - 1][i] + v[j + 1][i]);
                        v[j][i] = prev - w * (temp + f1(x[i], y[i])) / A;

                        //maxR1 += R[i][j] * R[i][j];
                        currentEps = Math.Abs(v[j][i] - prev);
                        if (currentEps > Eps_max)
                            Eps_max = currentEps;
                    }
                }

                p++;
                if ((Eps_max < Eps) || (p > N_max))
                    break;
            }
            // nevyazka na vyhode
            temp = 0.0;
            for (int j = 1; j < m; j++)
            {
                for (int i = 1; i < n; i++)
                {
                    temp = A * v[j][i] + h2 * (v[j][i - 1] + v[j][i + 1]) + k2 * (v[j - 1][i] + v[j + 1][i]) + f1(x[i], y[i]);
                    maxR1 += temp * temp;
                }
            }
            maxR1 = Math.Sqrt(maxR1);

            // table

            dataGridView1.Rows.Clear();
            dataGridView1.Columns.Clear();
            dataGridView1.Columns.Add("C1", "");
            dataGridView1.Columns[0].Width = 50;
            dataGridView1.Columns[0].Frozen = true;
            dataGridView1.Columns.Add("C2", "i");
            dataGridView1.Columns[1].Width = 50;
            dataGridView1.Columns[1].Frozen = true;

            dataGridView2.Rows.Clear();
            dataGridView2.Columns.Clear();
            dataGridView2.Columns.Add("C1", "");
            dataGridView2.Columns[0].Width = 50;
            dataGridView2.Columns[0].Frozen = true;
            dataGridView2.Columns.Add("C2", "i");
            dataGridView2.Columns[1].Width = 50;
            dataGridView2.Columns[1].Frozen = true;

            dataGridView3.Rows.Clear();
            dataGridView3.Columns.Clear();
            dataGridView3.Columns.Add("C1", "");
            dataGridView3.Columns[0].Width = 50;
            dataGridView3.Columns[0].Frozen = true;
            dataGridView3.Columns.Add("C2", "i");
            dataGridView3.Columns[1].Width = 50;
            dataGridView3.Columns[1].Frozen = true;

            for (int i = 0; i <= n; i++)                        //Создание столбцов для таблиц
            {
                dataGridView1.Columns.Add(Convert.ToString(buffer), Convert.ToString(buffer));
                dataGridView2.Columns.Add(Convert.ToString(buffer), Convert.ToString(buffer));
                dataGridView3.Columns.Add(Convert.ToString(buffer), Convert.ToString(buffer));

            }

            dataGridView1.Rows.Add("j", "Y\\X");  // Создание второй строки
            dataGridView2.Rows.Add("j", "Y\\X");  // Создание второй строки
            dataGridView3.Rows.Add("j", "Y\\X");  // Создание второй строки

            for (int i = 0; i <= n; i++)               //Заполнение второй строки
            {
                dataGridView1.Rows[0].Cells[i + 2].Value = x[i];//+2
                dataGridView2.Rows[0].Cells[i + 2].Value = x[i];
                dataGridView3.Rows[0].Cells[i + 2].Value = x[i];

            }
            for (int j = 0; j <= m; j++)          //Заполнение первых двух столбцов
            {
                dataGridView1.Rows.Add();
                dataGridView2.Rows.Add();
                dataGridView3.Rows.Add();

                for (int i = 0; i <= 1; i++)
                {
                    dataGridView1.Rows[j + 1].Cells[0].Value = j;
                    dataGridView1.Rows[j + 1].Cells[1].Value = y[j];
                    dataGridView2.Rows[j + 1].Cells[0].Value = j;
                    dataGridView2.Rows[j + 1].Cells[1].Value = y[j];
                    dataGridView3.Rows[j + 1].Cells[0].Value = j;
                    dataGridView3.Rows[j + 1].Cells[1].Value = y[j];
                }
            }
            double xMax = 0.0;
            double yMax = 0.0;
      
            for (int j = 0; j <= m; j++)              //Заполнение таблиц значениями
            {
                for (int i = 0; i <= n; i++)
                {
                    Pogr = Math.Abs(u[i][j] - v[i][j]);
                    v[i][j] = Math.Round(v[i][j] * 1000) / 1000;
                    u[i][j] = Math.Round(u[i][j] * 1000) / 1000;
                   
                    dataGridView1.Rows[j + 1].Cells[i + 2].Value = v[i][j];

                    dataGridView2.Rows[j + 1].Cells[i + 2].Value = u[i][j];
                   
                    dataGridView3.Rows[j + 1].Cells[i + 2].Value = Pogr;

                    if (Pogr > MaxPogr)
                    {
                        MaxPogr = Pogr;
                        xMax = x[i];
                        yMax = y[j];
                    }
                }
               
            }

            // Справка
            textBox9.Text = Convert.ToString(p);
            textBox10.Text = Convert.ToString(Eps_max);
            textBox11.Text = Convert.ToString(MaxPogr);
            textBox15.Text = Convert.ToString(MaxF);
            textBox16.Text = Convert.ToString(maxR1);
            textBox12.Text = Convert.ToString(xMax);
            textBox13.Text = Convert.ToString(yMax);

            textBox14.Text = "Нулевое начальноe приближение";

           
        }
    }
}
