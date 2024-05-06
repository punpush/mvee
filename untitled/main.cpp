#include <QGuiApplication>
#include <QQmlApplicationEngine>
#include <ellips.h>
#include <random>

int main(int argc, char *argv[])
{
    std::random_device dev;
    std::mt19937 rng(dev());
    std::normal_distribution<double> uint_dist10(0,10);

    std::vector<std::vector<double>> data {10, std::vector<double>(2, uint_dist10(rng))};
    Ellips el(data, 0.001);

    QGuiApplication app(argc, argv);

    QQmlApplicationEngine engine;
    const QUrl url(QStringLiteral("qrc:/untitled/Main.qml"));
    QObject::connect(
        &engine,
        &QQmlApplicationEngine::objectCreationFailed,
        &app,
        []() { QCoreApplication::exit(-1); },
        Qt::QueuedConnection);
    engine.load(url);

    return app.exec();
}
